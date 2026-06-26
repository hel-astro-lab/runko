# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Particle-in-cell simulation of driven kinetic turbulence using runko.
Driving is achived with oscillating Langevin antenna.

This setup assumes that:
- simulation box is a cube
- parallel and perpendiular components of each driven mode are the same
"""

import runko
import numpy as np
import itertools
import logging
import argparse

if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    logger = runko.runko_logger()

    if runko.on_main_rank():
        pass
        # logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description="PIC collisionless shock simulation")
    parser.add_argument("--conf", type=str, default="langevin_antenna_3d.ini",
                        help="Path to .ini configuration file")
    args = parser.parse_args()

    conf = runko.Configuration(args.conf)

    conf.n_laps = int(10 * conf.NxMesh * conf.Nx / conf.cfl)

    # --------------------------------------------------
    # Problem specific configuration

    Bcoeff = conf.Bperp_to_B0
    cold_sigma = conf.sigma

    n_e_num = conf.ppc
    modes_par, modes_perp = conf.modes_par, conf.modes_perp
    modes = [(modes_perp, 0, modes_par),
             (modes_perp, 0, -modes_par),
             (0, modes_perp, modes_par),
             (0, modes_perp, -modes_par)]
    delgam = conf.theta0 # temperature
    temp_ratio = conf.theta1_to_theta0 # T_i / T_e

    n_num = 2 * n_e_num

    plasma_freq_num = conf.cfl / conf.n_cells_per_skindepth
    average_gamma = 1 + 3 * delgam # approximation from Nättilä & Beloborodov (2021)
    conf.q0 = -average_gamma * (plasma_freq_num**2.0) / (n_e_num * (1.0 + conf.m0 / conf.m1))
    conf.q1 = abs(conf.q0)

    m0_num = conf.m0 * abs(conf.q0)
    m1_num = conf.m1 * abs(conf.q1)

    delgam0 = delgam
    delgam1 = temp_ratio * delgam0

    # Sigma that takes into account the heat contributions:
    sigma = cold_sigma / average_gamma
    B0_num = np.sqrt(n_num * (conf.cfl**2.0) * sigma * m0_num)

    # Assumes cubic box.
    L_num = conf.n_tiles[0] * conf.n_cells_per_tile[0]
    kpar_num = 2 * np.pi * modes_par / L_num
    kperp_num = 2 * np.pi * modes_perp / L_num

    A_num = kpar_num * B0_num / (kperp_num**2 * np.sqrt(len(modes)))

    valf = np.sqrt(sigma / (1 + sigma)) # in units of c
    linear_freq_num = kpar_num * conf.cfl * valf

    w0_num = 0.8 * linear_freq_num
    gamma0_num = -0.6 * linear_freq_num

    base_time_evolution = time_evolution = runko.sample_oscillating_langevin_antenna(size=conf.n_laps,
                                                                                     characteristic_freq=w0_num,
                                                                                     decorrelation_rate=gamma0_num,
                                                                                     gen=rng)


    def antenna(mode):
        # Note that in multi-node simulations, it is important to give this `gen` argument
        # with numpy rng that has a defined seed and every node uses the same seed.
        #
        # Otherwise, different parts of the modes will evolve differently.

        random_phase = np.exp(2j * np.pi * rng.random())
        if conf.fixed_antenna:
            return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=base_time_evolution * random_phase)
        else:
            time_evolution = runko.sample_oscillating_langevin_antenna(size=conf.n_laps,
                                                                       characteristic_freq=w0_num,
                                                                       decorrelation_rate=gamma0_num,
                                                                       gen=rng)
            return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=time_evolution * random_phase)


    antenna_modes = list(map(antenna, modes))
    # --------------------------------------------------
    # Print setup summary
    if runko.on_main_rank():
        from runko.auto_outdir import resolve_outdir
        W = 21
        logger.info(f"{'--- [io] ---':}")
        logger.info(f"  {'outdir':<{W}}= {conf.io_outdir}")
        logger.info(f"  {'resolved outdir':<{W}}= {resolve_outdir(conf)}")
        logger.info(f"  {'prefix / postfix':<{W}}= {conf.io_outidr_prefix} / {conf.io_outidr_postfix}")
        logger.info(f"  {'output_interval':<{W}}= {conf.io_output_interval}")

        logger.info(f"{'--- [grid] ---':}")
        logger.info(f"  {'tiles':<{W}}= {conf.n_tiles}")
        logger.info(f"  {'mesh per tile':<{W}}= {conf.n_cells_per_tile}")
        full_grid = np.array(conf.n_tiles) * np.array(conf.n_cells_per_tile)
        logger.info(f"  {'full grid':<{W}}= {full_grid}")

        logger.info(f"{'--- [simulation] ---':}")
        logger.info(f"  {'laps':<{W}}= {conf.n_laps}")
        logger.info(f"  {'cfl':<{W}}= {conf.cfl}")

        logger.info(f"{'--- [particles] ---':}")
        logger.info(f"  {'ppc':<{W}}= {n_e_num}")
        logger.info(f"  {'m0 / m1':<{W}}= {conf.m0} / {conf.m1}")
        logger.info(f"  {'q0 / q1':<{W}}= {conf.q0:.6g} / {conf.q1:.6g}")
        logger.info(f"  {'theta_e / theta_i':<{W}}= {delgam0:.6g} / {delgam1:.6g}")

        logger.info(f"{'--- [problem] ---':}")
        logger.info(f"  {'cold sigma':<{W}}= {cold_sigma}")
        logger.info(f"  {'sigma':<{W}}= {sigma}")
        logger.info(f"  {'B_init':<{W}}= {B0_num:.6g}")
        logger.info(f"  {'A':<{W}}= {A_num:.6g}")
        logger.info(f"  {'linear freq':<{W}}= {linear_freq_num:.6g}")
        logger.info(f"  {'n_filter_passes':<{W}}= {conf.n_filter_passes}")

        logger.info(f"{'--- [algorithms] ---':}")
        logger.info(f"  {'field_propagator':<{W}}= {conf.field_propagator}")
        logger.info(f"  {'particle_pusher':<{W}}= {conf.particle_pusher}")
        logger.info(f"  {'field_interpolator':<{W}}= {conf.field_interpolator}")
        logger.info(f"  {'current_depositer':<{W}}= {conf.current_depositer}")
        logger.info(f"  {'current_filter':<{W}}= {conf.current_filter}")
        logger.info(f"  {'tile_partitioning':<{W}}= {conf.tile_partitioning}")


    zero_field = lambda x, y, z: np.zeros_like(x)
    bz = lambda x, y, z: np.ones_like(x) * B0_num

    def pgen0(x, y, z):
        N = len(x)

        # Here we don't sample off rng as we want different mpi ranks to use different seeds.
        dx = np.random.random(N)
        dy = np.random.random(N)
        dz = np.random.random(N)

        # Particles 1 are going on top of particles 0,
        # so these positions has to be saved such that pgen1 can get them.
        pgen0.pos = x + dx, y + dy, z + dz
        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    tile_grid = runko.TileGrid(conf)

    if not tile_grid.initialized_from_restart_file():
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, conf)
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

    simulation = tile_grid.configure_simulation(conf)

    def sync_EB(x):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        x.comm_external(*EB)
        x.comm_local(*EB)

    simulation.prelude(sync_EB)

    # --------------------------------------------------
    # Main PIC loop

    def pic_simulation_step(x):

        # --- half B push + wall BC ---
        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- particle push + reflect + communicate ---
        x.prtcl_push()
        x.prtcl_pack_outgoing()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        if simulation.lap % 5 == 0:
            x.prtcl_sort()

        # --- current deposit + communicate ---
        x.prtcl_deposit_current()
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J)

        # --- current filter ---
        for i in range(conf.n_filter_passes):
            if i > 0 and i % 3 == 0:
                x.comm_external(runko.tools.comm_mode.emf_J)
                x.comm_local(runko.tools.comm_mode.emf_J)
            x.grid_filter_current()

        # --- second half B push + wall BC ---
        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        # --- E push + add current ---
        x.grid_push_e()
        x.grid_deposit_antenna_current()
        x.grid_add_current()
        x.comm_external(runko.tools.comm_mode.emf_E)
        x.comm_local(runko.tools.comm_mode.emf_E)

        # --- IO ---
        x.io_average_kinetic_energy()
        x.io_average_B_energy_density()
        x.io_average_E_energy_density()
        x.io_ram_usage()


        if simulation.lap % conf.io_output_interval == 0:
            x.io_emf_snapshot()
            x.io_prtcl_snapshot()
            simulation.log_timer_statistics()

    simulation.for_each_lap(pic_simulation_step)
    simulation.log_timer_statistics()
