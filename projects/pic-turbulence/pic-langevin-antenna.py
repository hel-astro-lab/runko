"""
Particle-in-cell simulation of driven kinetic turbulence using runko.
Driving is achived with oscillating Langevin antenna.
"""

import runko
import numpy as np
import itertools
import logging


if __name__ == "__main__":
    # Important for all of the nodes to have same seed.
    rng = np.random.default_rng(seed=42)
    logger = runko.runko_logger()
    config = runko.Configuration(None)

    # Input parameters:
    n_e_num = 8
    skin_depth_in_delta_x = 10
    modes_par, modes_perp = 1, 1
    modes = [(modes_perp, 0, modes_par),
             (modes_perp, 0, -modes_par),
             (0, modes_perp, modes_par),
             (0, modes_perp, -modes_par)]
    delgam = 0.3 # temperature
    temp_ration = 1 # T_i / T_e
    cold_sigma = 10 # magnetization (omega_ce/omega_pe)^2


    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 100
    config.NyMesh = 100
    config.NzMesh = 100
    config.cfl = 0.45
    config.Nt = int(10 * config.NxMesh * config.Nx / config.cfl)
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.field_propagator = "FDTD2"
    config.m0 = 1
    config.m1 = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"
    config.current_filter = "binomial2"
    config.tile_partitioning = "hilbert_curve"
    config.outdir = f"turb-langevin-antenna"


    # Problem specific configuration
    n_num = 2 * n_e_num

    plasma_freq_num = config.cfl / skin_depth_in_delta_x
    average_gamma = 1 + 3 * delgam # approximation from Nättilä & Beloborodov (2021)
    config.q0 = -average_gamma * (plasma_freq_num**2.0) / (n_e_num * (1.0 + config.m0 / config.m1))
    config.q1 = abs(config.q0)

    m0_num = config.m0 * abs(config.q0)
    m1_num = config.m1 * abs(config.q1)

    delgam0 = delgam
    delgam1 = temp_ration * delgam0

    # Sigma that takes into account the heat contributions:
    sigma = cold_sigma / average_gamma
    B0_num = np.sqrt(n_num * (config.cfl**2.0) * sigma * m0_num)

    Lx_num = config.Nx * config.NxMesh
    Ly_num = config.Ny * config.NyMesh
    kx_num = np.array(list(map(lambda n: 2 * np.pi * n[0] / Lx_num, modes)))
    ky_num = np.array(list(map(lambda n: 2 * np.pi * n[1] / Ly_num, modes)))
    logger.info(f"kx: {kx_num}")
    logger.info(f"ky: {ky_num}")
    logger.info(f"k_perp2: {np.sum(kx_num**2 + ky_num**2)}")

    A_num = np.sqrt(2 * B0_num**2 / np.sum(kx_num**2 + ky_num**2))

    Lz_num = config.Nz * config.NzMesh
    kpar_num = 2 * np.pi * modes_par / Lz_num

    # Linear frequency: kz * v_A
    linear_freq_num = kpar_num * config.cfl * np.sqrt(sigma / (1 + sigma))

    w0 = 0.8 * linear_freq_num
    gamma0 = -0.6 * linear_freq_num

    def antenna(mode):
        # Note that in multi-node simulations, it is important to give this `gen` argument
        # with numpy rng that has a defined seed and every node uses the same seed.
        #
        # Otherwise, different parts of the modes will evolve differently.
        time_evolution = runko.sample_oscillating_langevin_antenna(size=config.Nt,
                                                                   characteristic_freq=w0,
                                                                   decorrelation_rate=gamma0,
                                                                   gen=rng)
        random_phase = np.exp(2j * np.pi * rng.random())
        return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=time_evolution * random_phase)

    antenna_modes = list(map(antenna, modes))

    logger.info(f"Alfven vel: {np.sqrt(sigma / (1. + sigma))}")
    logger.info(f"cold sigma: {cold_sigma}")
    logger.info(f"sigma: {sigma}")
    logger.info(f"gamma_th: {average_gamma}")
    logger.info(f"B_guide: {B0_num}")
    logger.info(f"q0: {config.q0}")
    logger.info(f"q1: {config.q1}")
    logger.info(f"m0: {config.m0}")
    logger.info(f"m1: {config.m1}")
    logger.info(f"linear frequency: {linear_freq_num}")
    logger.info(f"characteristic freq: {w0}")
    logger.info(f"decorrelation rate: {gamma0}")
    logger.info(f"A0: {A_num}")

    zero_field = lambda x, y, z: np.zeros_like(x)
    bz = lambda x, y, z: np.ones_like(x) * B0_num

    def pgen0(x, y, z):
        N = len(x)

        dx = rng.random(N)
        dy = rng.random(N)
        dz = rng.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz
        # Here we don't sample off rng as we want different mpi ranks to use different seeds.a
        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    tile_grid = runko.TileGrid(config)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            tile.batch_set_EBJ(zero_field, zero_field, zero_field,
                               zero_field, zero_field, bz,
                               zero_field, zero_field, zero_field)
            for antenna in antenna_modes:
                logger.info("registering antenna...")
                tile.register_antenna(antenna)
            for _ in range(n_e_num):
                logger.info("Injecting particles of type 0...")
                tile.batch_inject_to_cells(0, pgen0)
                logger.info("Injecting particles of type 1...")
                tile.batch_inject_to_cells(1, pgen1)
            tile_grid.add_tile(tile, idx)

    # Initializes the simulation and returns a handle to it:
    simulation = tile_grid.configure_simulation(config)

    def sync_EB(x):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        x.comm_external(*EB)
        x.comm_local(*EB)

    simulation.prelude(sync_EB)

    def pic_simulation_step(x):

        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        x.prtcl_push()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            x.prtcl_sort()

        x.prtcl_deposit_current()
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)

        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J)

        x.grid_filter_current()
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J)
        x.grid_filter_current()
        x.grid_filter_current()

        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        x.grid_push_e()
        x.grid_deposit_antenna_current()
        x.grid_add_current()
        x.comm_external(runko.tools.comm_mode.emf_E)
        x.comm_local(runko.tools.comm_mode.emf_E)

        x.io_average_kinetic_energy()
        x.io_average_B_energy_density()
        x.io_average_E_energy_density()

        if simulation.lap % 20 == 0:
            x.io_emf_snapshot()

        if simulation.lap % 20 == 0:
            simulation.log_timer_statistics()

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
