import mpi_unittest
import itertools
import numpy as np

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
        particles.append(runko.ParticleState(pos=pos, vel=vel))

    return particles


def almost_equal(a: runko.ParticleState, b: runko.ParticleState, epsilon=0.00001):
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
        tile = runko.pic.Tile(idx, conf)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    vtiles = list(simulation.virtual_tiles())

    asserts = []
    for vtile in vtiles:
        asserts.append(mpi_unittest.assertEqualDeferred(type(vtile), runko.pic.Tile))

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
        tile = runko.pic.Tile(idx, conf)
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
                    s = runko.ParticleState(pos=(x, y, z), vel=(vx, vy, vz))
                    e = [p for p in expected_particles if almost_equal(p, s)]
                    asserts.append(mpi_unittest.assertEqualDeferred(1, len(e)))

        mpi_unittest.assertDeferredResults(asserts)
        mpi_unittest.assertEqual(True, len(asserts) >= 1)

    def f(local_tile, communicate, *_):
        communicate.virtual_tile_sync(runko.comm_mode.pic_particle)
        communicate.pairwise_moore(runko.comm_mode.pic_particle)

        local_tile.push_e()
        local_tile.push_half_b()

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
        tile = runko.pic.Tile(idx, conf)
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
                    s = runko.ParticleState(pos=(x, y, z), vel=(vx, vy, vz))
                    e = [p for p in expected_particles if almost_equal(p, s)]
                    asserts.append(mpi_unittest.assertEqualDeferred(1, len(e)))

        mpi_unittest.assertDeferredResults(asserts)
        mpi_unittest.assertEqual(True, len(asserts) >= 1)

    def f(local_tile, communicate, *_):
        communicate.virtual_tile_sync(runko.comm_mode.pic_particle)
        communicate.pairwise_moore(runko.comm_mode.pic_particle)

        local_tile.push_e()
        local_tile.push_half_b()

    simulation.for_one_lap(f)
    assertChangedParticles()
    simulation.for_each_lap(f)
    assertChangedParticles()


if __name__ == "__main__":
    virtual_pic_tiles()
    pic_noop_communication()
    pic_communication()
