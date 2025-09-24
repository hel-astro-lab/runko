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

    config.cfl = 0.45
    config.Nt = int(1000 / config.cfl)
    config.Nx = 8
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 50
    config.NyMesh = 50
    config.NzMesh = 250
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.field_propagator = "FDTD2"
    config.m0 = 1
    config.m1 = config.m0
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"
    config.current_filter = "binomial2"
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    if len(sys.argv) > 1:
        config.outdir = sys.argv[1]

    # Start of problem specific configuration

    # Free variables:
    # Plasma motion in lab frame
    Gamma = 10
    # Plasma magnetisation in co-moving fluid rest frame
    sigma = 1
    # Cells-per-skindepth perpendicular to bulk motion direction in both rest and lab frames (c_omp = C / OMega Plasma)
    c_omp = 5
    # Plasma density per species
    ppc = 32
    # Plasma temperature in co-moving fluid rest frame in units of electron rest mass
    temperature = 0.3

    # Dependent variables:
    # Plasma motion (beta gamma factor)
    v_x = -np.sqrt(1-1/Gamma**2)
    u_x = Gamma*v_x
    # Plasma density (overall, all species per cell)
    oppc = 2*ppc
    # Particle charge
    config.q0 = -config.cfl * c_omp * np.sqrt(config.m0 * Gamma / oppc)
    config.q1 = -config.q0
    # Magnetic field strength in the lab frame
    B_z = config.cfl * np.sqrt(config.m0 * Gamma * sigma * oppc)
    # Electric field strength in the lab frame
    E_y = B_z

    # End of problem specific configuration

    logger.info(f"Bulk Gamma and v_x: {Gamma}, {v_x}")
    logger.info(f"Particle thermal spread: {delgam}")
    logger.info(f"Plasma beta: absent")
    logger.info(f"Plasma sigma: {sigma}")
    logger.info(f"B_z, E_y: {B_z}, {E_y}")
    logger.info(f"q0: {config.q0}")
    logger.info(f"q1: {config.q1}")
    logger.info(f"m0: {config.m0}")
    logger.info(f"m1: {config.m1}")
    
    # Configure fields
    zero_field = lambda x, y, z: np.zeros_like(x)

    Bx = zero_field
    By = zero_field
    Bz = lambda x, y, z: np.full_like(x, B_z)

    Ex = zero_field
    Ey = lambda x, y, z: np.full_like(x, E_y)
    Ez = zero_field

    # Configure particles
    def pgen0(x, y, z):
        N = len(x)

        dx = rng.random(N)
        dy = rng.random(N)
        dz = rng.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz
        vel = runko.sample_boosted_juttner_synge(N, temperature, gen=rng, Gamma=Gamma, direction='-x')
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), temperature, gen=rng, Gamma=Gamma, direction='-x')
        return runko.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    # Initialise tiles

    # TileGrid ctor:
    # - balances tiles based on conf (catepillar, hilbert)
    # - checks if restarts files are present for current config
    #   - if found initializes the tiles based on the restart files
    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.Tile(idx, config)
            tile.batch_set_EBJ(Ex, Ey, Ez,
                               Bx, By, Bz,
                               zero_field, zero_field, zero_field)
            for _ in range(ppc):
                runko.runko_logger().info("Injecting particles of type 0...")
                tile.batch_inject_to_cells(0, pgen0)
                runko.runko_logger().info("Injecting particles of type 1...")
                tile.batch_inject_to_cells(1, pgen1)
            tile_grid.add_tile(tile, idx)

    # Initialises the simulation and returns a handle to it:
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
