"""
Particle-in-cell simulation of driven kinetic turbulence using runko.
Driving is achived with oscillating Langevin antenna.
"""

import runko
import numpy as np
import itertools
import logging


if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    logger = runko.runko_logger()

    if runko.on_main_rank():
        pass
        # logger.setLevel(logging.DEBUG)

    config = runko.Configuration(None)

    config.outdir = "turb-langevin-antenna"
    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 100
    config.NyMesh = 100
    config.NzMesh = 100
    config.cfl = 0.45
    config.Nt = int(5.0*config.Nx*config.NxMesh/config.cfl/3.0) # =5eddy turnover times
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
    sigma = 10 # magnetization (omega_ce/omega_pe)^2

    delgam0 = delgam
    delgam1 = temp_ration * delgam0

    # No corrections; cold sigma
    binit_nc = np.sqrt(oppc * (config.cfl**2.0) * sigma * m0)
    # another approximation which is more accurate at \delta ~ 1
    gammath = 1.0 + (3.0/2.0) * delgam1

    binit_approx = np.sqrt(gammath * oppc * m0 * (config.cfl**2.0) * sigma)
    binit = binit_approx # NOTE: selecting this as our sigma definitions

    # Antenna setup (Tenbarge et. al. 2014):
    drive_amplitude = 0.8 * binit
    drive_characteristic_freq = 0.8 # in units of linear frequency
    drive_decorrelation_rate = 0.8 # in units of linear frequency

    mode_number = 1
    eddy_length = config.Nz * config.NzMesh / mode_number
    light_crossing_time = eddy_length / config.cfl
    v_alfven = config.cfl * np.sqrt(sigma / (1.0 + sigma))
    kz_alfven = 2 * np.pi / eddy_length
    linear_freq = kz_alfven * v_alfven

    w0 = drive_characteristic_freq * linear_freq
    gamma0 = -drive_decorrelation_rate * linear_freq

    modes = (lambda N: [(N, 0, N), (N, 0, -N), (0, N, N), (0, N, -N)])(mode_number)
    k_perp = 2 * np.pi * mode_number / (config.Nx * config.NxMesh)
    A0 = config.cfl * kz_alfven * binit * np.sqrt(2) / (v_alfven * k_perp**2 * np.sqrt(len(modes)))

    # Note that in multi-node simulations, # it is important to give this `gen` argument
    # with numpy rng that has a defined seed # and every node uses the same seed.
    #
    # Otherwise, different parts of the modes will evolve differently.
    # Here it does not matter, because we only have one node.
    base_time_evolution = runko.sample_oscillating_langevin_antenna(size=config.Nt,
                                                                    characteristic_freq=w0,
                                                                    decorrelation_rate=gamma0,
                                                                    gen=rng)

    random_phases = [np.exp(2j * np.pi * rng.random()) for _ in modes]
    time_evolutions = [phi * base_time_evolution for phi in random_phases]
    antenna_mode = runko.emf.threeD.antenna_mode
    antenna_modes = [antenna_mode(A=(0, 0, A0), n=mode, lap_coeffs=time) for mode, time in zip(modes, time_evolutions)]

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
    logger.info(f"linear frequency: {linear_freq}")
    logger.info(f"characteristic freq: {w0}")
    logger.info(f"decorrelation rate: {gamma0}")
    logger.info(f"A0: {A0}")

    zero_field = lambda x, y, z: np.zeros_like(x)
    bz = lambda x, y, z: np.ones_like(x) * binit

    def pgen0(x, y, z):
        N = len(x)

        dx = rng.random(N)
        dy = rng.random(N)
        dz = rng.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz
        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0, gen=rng)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0, gen=rng)
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
            for _ in range(ppc):
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
