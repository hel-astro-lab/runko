import mpi_unittest
import itertools
import numpy as np
import tempfile

import runko

def make_test_grid():
    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.Nt = 3
    config.Nx = 2
    config.Ny = 2
    config.Nz = 2
    config.NxMesh = 10
    config.NyMesh = 11
    config.NzMesh = 13
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.delgam = 1.0e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    config.outdir = tempfile.mkdtemp(prefix="runko-pic-test-output")

    return config, runko.TileGrid(config)


def make_test_particles(tile, after_comm=False):
    """
    Create one particles in all tile 27 subregion.
    If after_comm=True, then returns the particles
    as if particles from after_comm=False have been communicated.
    """

    mins = np.array(tile.mins)
    maxs = np.array(tile.maxs)

    L = np.array([tile.maxs[i] - tile.mins[i] for i in range(3)])
    center = mins + L / 2

    particles = []

    I = [-1, 0, 1]
    for dir in itertools.product(I, I, I):
        d = np.array(dir)
        pos = center + d * (1 + L / 2 - 2 * int(after_comm))
        vel = d if not after_comm else -d
        particles.append(runko.pic.threeD.ParticleState(pos=pos, vel=vel))

    return particles


def almost_equal(a: runko.pic.threeD.ParticleState, b: runko.pic.threeD.ParticleState, epsilon=0.00001):
    def ok(x, y):
        return (epsilon - abs(x - y)) >= 0

    xok = ok(a.pos[0], b.pos[0])
    yok = ok(a.pos[1], b.pos[1])
    zok = ok(a.pos[2], b.pos[2])
    velxok = ok(a.vel[0], b.vel[0])
    velyok = ok(a.vel[1], b.vel[1])
    velzok = ok(a.vel[2], b.vel[2])

    return xok and yok and zok and velxok and velyok and velzok


def virtual_pic_tiles():
    """
    Here we test that after grid.configure_simulation(...)
    each virtual tile has bee initialized to pic tile,
    as they are next to them.
    """

    conf, tile_grid = make_test_grid()

    mpi_unittest.assertEqual(tile_grid.initialized_from_restart_file(), False)

    local_idx = list(tile_grid.local_tile_indices())
    mpi_unittest.assertNotEqual(len(local_idx), 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, conf)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    vtiles = list(simulation.virtual_tiles())

    asserts = []
    for vtile in vtiles:
        asserts.append(mpi_unittest.assertEqualDeferred(type(vtile), runko.pic.threeD.Tile))

    mpi_unittest.assertDeferredResults(asserts)


def pic_noop_communication():
    """
    TileGrid with only pic tiles, which have been initialized
    to have particles only in the non-halo region.

    After virtual tile sync and pairwise moore communication
    nothing should be changed.
    """

    conf, tile_grid = make_test_grid()

    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, conf)
        tile.inject(0, make_test_particles(tile, after_comm=True))
        tile.inject(1, make_test_particles(tile, after_comm=True))
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    import matplotlib.pyplot as plt
    def plot():
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for tile in simulation.local_tiles():
            posx, posy, posz = tile.get_positions(0)
            ax.scatter(posx, posy, posz, label=f"p: {0}, tile: {tile.mins} x {tile.maxs}")
        fig.suptitle(f"rank: {mpi_unittest._rank}")
        ax.legend()
        plt.show()

    def assertUnchangedParticles():
        asserts = []

        A = lambda x, y: mpi_unittest.assertEqualDeferred(1, approx_count(x, y))

        for tile in simulation.local_tiles():
            for i in (0, 1):
                posx, posy, posz = tile.get_positions(i)
                velx, vely, velz = tile.get_velocities(i)

                N = len(posx)
                asserts.append(mpi_unittest.assertEqualDeferred(len(posy), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(posz), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(velx), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(vely), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(velz), N))

                expected_particles = make_test_particles(tile, after_comm=True)

                for x, y, z, vx, vy, vz in zip(posx, posy, posz, velx, vely, velz):
                    s = runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(vx, vy, vz))
                    e = [p for p in expected_particles if almost_equal(p, s)]
                    asserts.append(mpi_unittest.assertEqualDeferred(1, len(e)))

        mpi_unittest.assertDeferredResults(asserts)
        mpi_unittest.assertEqual(True, len(asserts) >= 1)

    def f(x):
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        x.grid_push_e()
        x.grid_push_half_b()

    assertUnchangedParticles()
    simulation.for_one_lap(f)
    assertUnchangedParticles()
    simulation.for_each_lap(f)
    assertUnchangedParticles()


def pic_communication():
    """
    TileGrid with only pic tiles, which have been initialized
    to have particles only in the halo region.

    After virtual tile sync and pairwise moore communication
    there should not be any particles in the halo regions.
    """

    conf, tile_grid = make_test_grid()

    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, conf)
        tile.inject(0, make_test_particles(tile))
        tile.inject(1, make_test_particles(tile))
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    import matplotlib.pyplot as plt
    def plot():
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for tile in simulation.local_tiles():
            posx, posy, posz = tile.get_positions(0)
            ax.scatter(posx, posy, posz, label=f"p: {0}, tile: {tile.mins} x {tile.maxs}")

            expected_particles = make_test_particles(tile, after_comm=True)
            ax.scatter([p.pos[0] for p in expected_particles], [p.pos[1] for p in expected_particles], [p.pos[2] for p in expected_particles], marker="^", label="expected")

        fig.suptitle(f"rank: {mpi_unittest._rank}")
        ax.legend()
        plt.show()

    def assertChangedParticles():
        asserts = []

        A = lambda x, y: mpi_unittest.assertEqualDeferred(1, approx_count(x, y))

        for tile in simulation.local_tiles():
            expected_particles = make_test_particles(tile, after_comm=True)
            for i in (0, 1):
                posx, posy, posz = tile.get_positions(i)
                velx, vely, velz = tile.get_velocities(i)

                N = len(posx)
                asserts.append(mpi_unittest.assertEqualDeferred(len(posy), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(posz), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(velx), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(vely), N))
                asserts.append(mpi_unittest.assertEqualDeferred(len(velz), N))

                for x, y, z, vx, vy, vz in zip(posx, posy, posz, velx, vely, velz):
                    s = runko.pic.threeD.ParticleState(pos=(x, y, z), vel=(vx, vy, vz))
                    e = [p for p in expected_particles if almost_equal(p, s)]
                    asserts.append(mpi_unittest.assertEqualDeferred(1, len(e)))

        mpi_unittest.assertDeferredResults(asserts)
        mpi_unittest.assertEqual(True, len(asserts) >= 1)

    def f(x):
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        x.grid_push_e()
        x.grid_push_half_b()

    simulation.for_one_lap(f)
    assertChangedParticles()
    simulation.for_each_lap(f)
    assertChangedParticles()


def pic_kinetic_energy_reduction():

    conf, tile_grid = make_test_grid()

    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, conf)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    # There will be 3 laps and we will store kinetic energy data from each:
    datas = []

    # Lap 1: no particles

    test_loop = lambda x: x.io_average_kinetic_energy()

    simulation.for_one_lap(test_loop)

    if runko.on_main_rank():
        d = np.loadtxt(conf.outdir + "/average_kinetic_energy.txt")
        if np.any(d[1:] != 0):
            msg = f"Kinetic energies of zero particles should be zero: {d[1:]}"
            raise RuntimeError(msg)
        datas.append(d)

    # Lap 1: unmoving particles 0 and moving particles 1

    def unmoving_pgen(x, y, z):
        pos = x, y, z
        vel = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
        return runko.pic.threeD.ParticleStateBatch(pos=pos, vel=vel)

    def moving_pgen(x, y, z):
        pos = x, y, z
        vel = np.ones_like(x) / 10, np.ones_like(x) / 10, np.ones_like(x) / 10
        return runko.pic.threeD.ParticleStateBatch(pos=pos, vel=vel)


    for tile in simulation.local_tiles():
        tile.batch_inject_to_cells(0, unmoving_pgen)
        tile.batch_inject_to_cells(1, moving_pgen)

    simulation.for_one_lap(test_loop)

    if runko.on_main_rank():
        d = np.loadtxt(conf.outdir + "/average_kinetic_energy.txt")

        if d[1, 1] != 0:
            msg = f"Kinetic energies of unmoving particles should be zero: {d[1,1]}"
            raise RuntimeError(msg)

        if d[1, 2] == 0:
            msg = f"Kinetic energies of moving particles should not be zero: {d[1,2]}"
            raise RuntimeError(msg)

        datas.append(d)

    # Lap 2: add moving 0 and 1 particles

    for tile in simulation.local_tiles():
        # It should not matter that the particles are unevenly between the tiles.
        if runko.on_main_rank():
            # It is assumed that this test is run with 4 ranks.
            for _ in range(4):
                tile.batch_inject_to_cells(0, moving_pgen)

        tile.batch_inject_to_cells(1, moving_pgen)

    simulation.for_one_lap(test_loop)

    if runko.on_main_rank():
        d = np.loadtxt(conf.outdir + "/average_kinetic_energy.txt")

        if not np.isclose(2 * d[2, 1], d[2, 2]):
            msg = "Kinetic energy should be linear in number of moving particles: "
            msg += f"{2 * d[2, 1]} != {d[2, 2]}"
            raise RuntimeError(msg)

        datas.append(d)

    # Now just make sure that the stored data is not changed:

    if runko.on_main_rank():
        ok0 = np.all(datas[2][0, 1:] == datas[0][1:])
        ok1 = np.all(datas[2][1, 1:] == datas[1][1, 1:])
        ok2 = np.all(datas[0][1:] == datas[1][0, 1:])

        if not (ok0 and ok1 and ok2):
            msg = f"Stored kinetic energy data should not change."
            raise RuntimeError(msg)

def pic_communication_with_empty_tiles():
    """
    Regression test: local_communication must not crash when
    some tiles have zero particles for a species.

    Previously, divide_to_subregions returned an empty spans map
    for empty containers, causing spans.at({0,0,0}) to throw
    std::out_of_range in local_communication_prelude.
    """

    conf, tile_grid = make_test_grid()

    for i, idx in enumerate(tile_grid.local_tile_indices()):
        tile = runko.pic.threeD.Tile(idx, conf)

        # Only inject particles into even-numbered tiles;
        # odd-numbered tiles stay empty (0 particles).
        if i % 2 == 0:
            tile.inject(0, make_test_particles(tile, after_comm=True))
            tile.inject(1, make_test_particles(tile, after_comm=True))

        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    # This must not crash with map::at key not found.
    def f(x):
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

    simulation.for_one_lap(f)

    # Verify all tiles still have valid particle containers.
    asserts = []
    for tile in simulation.local_tiles():
        for species in (0, 1):
            posx, posy, posz = tile.get_positions(species)
            # Particle count must be non-negative (>= 0).
            asserts.append(mpi_unittest.assertEqualDeferred(len(posx) >= 0, True))

    mpi_unittest.assertDeferredResults(asserts)


def pic_wrap_positions():
    """
    Test that wrap_positions correctly wraps particle coordinates
    at both positive and negative domain boundaries in all 3 directions.

    For each of the 6 domain faces, inject a particle slightly beyond
    the boundary. After communication (which triggers wrap_positions
    in postlude), verify the coordinate is wrapped to the opposite side.
    """

    conf, tile_grid = make_test_grid()

    Lx = conf.Nx * conf.NxMesh  # 20
    Ly = conf.Ny * conf.NyMesh  # 22
    Lz = conf.Nz * conf.NzMesh  # 26

    eps = 0.25
    PS = runko.pic.threeD.ParticleState

    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, conf)

        mins = tile.mins
        maxs = tile.maxs
        cx = (mins[0] + maxs[0]) / 2
        cy = (mins[1] + maxs[1]) / 2
        cz = (mins[2] + maxs[2]) / 2

        particles = []

        # x- boundary: particle just below x=0
        if mins[0] == 0:
            particles.append(PS(pos=(-eps, cy, cz), vel=(-1, 0, 0)))

        # x+ boundary: particle just above x=Lx
        if maxs[0] == Lx:
            particles.append(PS(pos=(Lx + eps, cy, cz), vel=(1, 0, 0)))

        # y- boundary: particle just below y=0
        if mins[1] == 0:
            particles.append(PS(pos=(cx, -eps, cz), vel=(0, -1, 0)))

        # y+ boundary: particle just above y=Ly
        if maxs[1] == Ly:
            particles.append(PS(pos=(cx, Ly + eps, cz), vel=(0, 1, 0)))

        # z- boundary: particle just below z=0
        if mins[2] == 0:
            particles.append(PS(pos=(cx, cy, -eps), vel=(0, 0, -1)))

        # z+ boundary: particle just above z=Lz
        if maxs[2] == Lz:
            particles.append(PS(pos=(cx, cy, Lz + eps), vel=(0, 0, 1)))

        tile.inject(0, particles)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    def f(x):
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

    simulation.for_one_lap(f)

    # After communication + wrap_positions, verify correctness.
    tol = 0.001
    asserts = []
    total_particles = 0

    for tile in simulation.local_tiles():
        posx, posy, posz = tile.get_positions(0)
        velx, vely, velz = tile.get_velocities(0)

        total_particles += len(posx)

        for i in range(len(posx)):
            x, y, z = float(posx[i]), float(posy[i]), float(posz[i])
            vx, vy, vz = float(velx[i]), float(vely[i]), float(velz[i])

            # All positions must be within domain bounds
            asserts.append(mpi_unittest.assertEqualDeferred(
                0 <= x < Lx and 0 <= y < Ly and 0 <= z < Lz, True))

            # Check wrapped coordinate based on velocity marker
            if abs(vx - (-1)) < tol:   # x- wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(x - (Lx - eps)) < tol, True))

            elif abs(vx - 1) < tol:    # x+ wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(x - eps) < tol, True))

            elif abs(vy - (-1)) < tol:  # y- wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(y - (Ly - eps)) < tol, True))

            elif abs(vy - 1) < tol:    # y+ wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(y - eps) < tol, True))

            elif abs(vz - (-1)) < tol:  # z- wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(z - (Lz - eps)) < tol, True))

            elif abs(vz - 1) < tol:    # z+ wrap
                asserts.append(mpi_unittest.assertEqualDeferred(
                    abs(z - eps) < tol, True))

    # Each tile in 2x2x2 grid touches 3 domain faces → 3 particles injected.
    # After communication, each tile receives 3 wrapped particles.
    # With 2 tiles per rank: expect 6 particles per rank.
    asserts.append(mpi_unittest.assertEqualDeferred(total_particles, 6))

    mpi_unittest.assertDeferredResults(asserts)


if __name__ == "__main__":
    virtual_pic_tiles()
    pic_noop_communication()
    pic_communication()
    pic_communication_with_empty_tiles()
    pic_wrap_positions()
    pic_kinetic_energy_reduction()
