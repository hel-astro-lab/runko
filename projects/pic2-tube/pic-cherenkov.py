"""
Runko particle-in-cell simulation studing the numerical Cherenkov instability.
"""

import matplotlib # This prevents a segfault when importing runko for some reason
import runko
import numpy as np
import itertools
import logging
import sys


if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    logger = runko.runko_logger()

    if runko.on_main_rank():
        pass
        # logger.setLevel(logging.DEBUG)

    config = runko.Configuration(None)

    config.Nt = 200
    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 80
    config.NyMesh = 80
    config.NzMesh = 80
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 0.45
    config.field_propagator = "FDTD2"
    config.m0 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"
    config.current_filter = "binomial2"
    config.tile_partitioning = "hilbert_curve"
    if len(sys.argv) > 1:
        config.outdir = sys.argv[1]

    # Problem specific configuration
    ppc = 2 # particles per cell (one particle type)
    oppc = 2 * ppc # overall particles per cell (all particle types)
    gamma = 1
    c_omp = 1
    omp = config.cfl / c_omp

    config.q0 = -gamma * (omp**2.0) / (0.5 * oppc * (1.0 + config.m0 / config.m1))
    config.q1 = abs(config.q0)

    m0 = config.m0 * abs(config.q0)
    m1 = config.m1 * abs(config.q1)

    delgam = 0.3 # temperature
    temp_ration = 1 # T_i / T_e
    sigma = 10 # temperature corrections
               # or magnetization number (omega_ce/omega_pe)^2, including gamma for inertia??

    delgam0 = delgam
    delgam1 = temp_ration * delgam0

    # No corrections; cold sigma
    binit_nc = np.sqrt(oppc * (config.cfl**2.0) * sigma * m0)
    # another approximation which is more accurate at \delta ~ 1
    gammath = 1.0 + (3.0/2.0) * delgam1

    binit_approx = np.sqrt(gammath * oppc * m0 * (config.cfl**2.0) * sigma)
    binit = binit_approx # NOTE: selecting this as our sigma definitions

    logger.info(f"Positron thermal spread: {delgam1}")
    logger.info(f"Electron thermal spread: {delgam0}")
    logger.info(f"Alfven vel: {np.sqrt(sigma / (1. + sigma))}")
    ion_beta = 2. * delgam1 / (sigma * (m1 /m0 + 1.) / (m1 /m0 ))
    logger.info(f"Ion beta: {ion_beta}")
    logger.info(f"Electron beta: {2.*delgam0/(sigma*(1/m0+1.))}")
    logger.info(f"sigma: {sigma}")
    logger.info(f"mass term: {np.sqrt(m0 + m1)}")
    logger.info(f"gamma_th: {gammath}")
    logger.info(f"B_guide (no corr): {binit_nc}")
    logger.info(f"B_guide (approx): {binit_approx}")
    logger.info(f"B_guide (used): {binit}")
    logger.info(f"q0: {config.q0}")
    logger.info(f"q1: {config.q1}")
    logger.info(f"m0: {config.m0}")
    logger.info(f"m1: {config.m1}")

    # Decaying setup:

    A0 = 0.8 * binit
    n_perp = 1
    n_par = 2

    # These use legacy rand for parity with pic-turbulence.
    np.random.seed(n_perp)
    ph1 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)
    ph2 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)
    ph3 = 2.0 * np.pi * np.random.rand(n_perp, n_perp, n_par)

    Lx = config.Nx * config.NxMesh / c_omp
    Ly = config.Ny * config.NyMesh / c_omp
    Lz = config.Nz * config.NzMesh / c_omp

    kx = 2.0 * np.pi / (c_omp * Lx)
    ky = 2.0 * np.pi / (c_omp * Ly)
    kz = 2.0 * np.pi / (c_omp * Lz)

    # Length of the wave mode vector
    def beta(n, m):
        return np.sqrt(8.0) / (np.sqrt(n**2 + m**2) * n_perp * np.sqrt(n_par))

    def Bx(x, y, z):
        bx = np.zeros_like(x)
        x, y = y, x # done indirectly in pic-trubulence/antenna3d.py
        I = range(1, n_perp + 1)
        for n, m in itertools.product(I, I):
            norm = beta(n ,m)
            for o in range(1, n_par + 1):
                xmodx = np.sin(m * kx * x + ph1[n - 1, m - 1, o - 1])
                xmody = np.cos(n * ky * y + ph2[n - 1, m - 1, o - 1])
                zmod = np.sin(o * kz * z + ph3[n - 1, m - 1, o - 1])
                bx += norm * n * xmodx * xmody * zmod

        return A0 * bx

    def By(x, y, z):
        by = np.zeros_like(x)
        x, y = y, x # done indirectly in pic-trubulence/antenna3d.py
        I = range(1, n_perp + 1)
        for n, m in itertools.product(I, I):
            norm = beta(n ,m)
            for o in range(1, n_par + 1):
                ymodx = np.cos(m * kx * x + ph1[n - 1, m - 1, o - 1])
                ymody = np.sin(n * ky * y + ph2[n - 1, m - 1, o - 1])
                zmod = np.sin(o * kz * z + ph3[n - 1, m - 1, o - 1])
                by -= norm * m * ymodx * ymody * zmod

        return A0 * by

    Bz = lambda x, y, z: np.full_like(x, binit)
    zero_field = lambda x, y, z: np.zeros_like(x)


    def pgen0(x, y, z):
        N = len(x)

        dx = rng.random(N)
        dy = rng.random(N)
        dz = rng.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz
        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0, gen=rng)
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0, gen=rng)
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    # TileGrid ctor:
    # - balances tiles based on conf (catepillar, hilbert)
    # - checks if restarts files are present for current config
    #   - if found initializes the tiles based on the restart files
    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.Tile(idx, config)
            tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                               Bx, By, Bz,
                               zero_field, zero_field, zero_field)
            for _ in range(ppc):
                runko.runko_logger().info("Injecting particles of type 0...")
                tile.batch_inject_to_cells(0, pgen0)
                runko.runko_logger().info("Injecting particles of type 1...")
                tile.batch_inject_to_cells(1, pgen1)
            tile_grid.add_tile(tile, idx)

    # Initializes the simulation and returns a handle to it:
    simulation = tile_grid.configure_simulation(config)

    def sync_EB(tile, comm, io):
        EB = (runko.comm_mode.emf_E, runko.comm_mode.emf_B)
        comm.virtual_tile_sync(*EB)
        comm.pairwise_moore(*EB)

    simulation.prelude(sync_EB)

    def pic_simulation_step(tile, comm, io):

        tile.push_half_b()
        comm.virtual_tile_sync(runko.comm_mode.emf_B)
        comm.pairwise_moore(runko.comm_mode.emf_B)

        tile.push_particles()
        comm.virtual_tile_sync(runko.comm_mode.pic_particle)
        comm.pairwise_moore(runko.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            tile.sort_particles()

        tile.deposit_current()
        comm.virtual_tile_sync(runko.comm_mode.emf_J)
        comm.pairwise_moore(runko.comm_mode.emf_J_exchange)

        comm.virtual_tile_sync(runko.comm_mode.emf_J)
        comm.pairwise_moore(runko.comm_mode.emf_J)
        tile.filter_current()
        comm.virtual_tile_sync(runko.comm_mode.emf_J)
        comm.pairwise_moore(runko.comm_mode.emf_J)
        tile.filter_current()
        tile.filter_current()

        tile.push_half_b()
        comm.virtual_tile_sync(runko.comm_mode.emf_B)
        comm.pairwise_moore(runko.comm_mode.emf_B)

        tile.push_e()
        tile.subtract_J_from_E()
        comm.virtual_tile_sync(runko.comm_mode.emf_E)
        comm.pairwise_moore(runko.comm_mode.emf_E)

        if simulation.lap % 20 == 0:
            io.emf_snapshot()

        if simulation.lap % 10 == 0:
            simulation.log_timer_statistics()


    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
