# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Particle-in-cell simulation of kinetic turbulence using runko.

Energy enters the box in one of two ways, selected with `use_antenna` in the .ini:
- True:  driven turbulence. B starts as a uniform guide field and an oscillating Langevin
         antenna injects energy at the driving scale.
- False: decaying turbulence. A superposition of sheared-Alfvén-like modes
         is stamped into Bx/By at t=0 on top of the guide field, and decays.
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

    parser = argparse.ArgumentParser(description="PIC turbulence simulation")
    parser.add_argument("--conf", type=str, default="turb_driven.ini",
                        help="Path to .ini configuration file")
    args = parser.parse_args()
    conf = runko.Configuration(args.conf)

    # --------------------------------------------------
    # Problem specific configuration

    L_num = conf.n_tiles[0] * conf.n_cells_per_tile[0] # assumes cubic box
    modes_par, modes_perp = conf.modes_par, conf.modes_perp

    n_e_num = conf.ppc # electron number density (ppc)
    n_num = 2 * n_e_num # total particle density 
    plasma_freq_num = conf.cfl / conf.n_cells_per_skindepth
    conf.q0 = -(plasma_freq_num**2.0) / n_num 
    conf.q1 = abs(conf.q0)
    m0_num = conf.m0 * abs(conf.q0)
    m1_num = conf.m1 * abs(conf.q1)

    # plasma enthalphy per particle in units of m c^2 (dimensionless);
    # fit formula for K_3(1/theta)/K_2(1/theta)
    delgam0 = conf.theta0
    delgam1 = conf.theta1_to_theta0 * delgam0
    enth_e = (1 + 3.89*delgam0 + 5.56*delgam0**2)/(1 + 1.39*delgam0)
    enth_i = (1 + 3.89*delgam1 + 5.56*delgam1**2)/(1 + 1.39*delgam1)

    sigma_e = conf.sigma # cold electron (species 0) sigma
    sigma_i = sigma_e * (m0_num/m1_num) # cold ion sigma (species 1)

    # magnetic field normalization based on cold electron sigma 
    B0_num = np.sqrt(n_e_num * (conf.cfl**2.0) * sigma_e * m0_num)

    # magnetization accounting for enthalphy
    sigma_e_hot = sigma_e/enth_e
    sigma_i_hot = sigma_i/enth_i

    # total hot sigma
    sigma = B0_num**2/(n_e_num*( enth_e*m0_num + enth_i*m1_num )*conf.cfl**2)

    # Alfven speed
    valf = np.sqrt(sigma / (1 + sigma)) # in units of c

    # driving-scale eddy is one perpendicular wavelength across the box
    # turnover is its light crossing time  
    eddy_length_num = L_num / modes_perp
    laps_per_eddy = eddy_length_num / conf.cfl
    conf.n_laps = int(conf.n_eddy_turnovers * laps_per_eddy)

    zero_field = lambda x, y, z: np.zeros_like(x)
    bz = lambda x, y, z: np.ones_like(x) * B0_num

    if conf.use_antenna: # driven: energy injection via oscillating Langevin antenna
        modes = [(modes_perp, 0, modes_par),
                 (modes_perp, 0, -modes_par),
                 (0, modes_perp, modes_par),
                 (0, modes_perp, -modes_par)]

        kpar_num = 2 * np.pi * modes_par / L_num
        kperp_num = 2 * np.pi * modes_perp / L_num

        A_num = kpar_num * B0_num*conf.Bperp_to_B0 / (kperp_num**2 * np.sqrt(len(modes)))

        linear_freq_num = kpar_num * conf.cfl * valf

        w0_num = conf.antenna_freq * linear_freq_num # omega_0 for antenna
        gamma0_num = -conf.antenna_decorr * linear_freq_num # omega_dec for antenna

        def antenna_mode_for(mode):
            # Note that in multi-node simulations, it is important to give this `gen` argument
            # with numpy rng that has a defined seed and every node uses the same seed.
            # Otherwise, different parts of the modes will evolve differently.

            random_phase = np.exp(2j * np.pi * rng.random())
            time_evolution = runko.sample_oscillating_langevin_antenna(size=conf.n_laps,
                                                                       characteristic_freq=w0_num,
                                                                       decorrelation_rate=gamma0_num,
                                                                       gen=rng)
            return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=time_evolution * random_phase)

        antenna_modes = list(map(antenna_mode_for, modes))
        bx = by = zero_field
    else: # decaying: all energy is put into the initial perpendicular field
        antenna_modes = []

        A0_num = B0_num*conf.Bperp_to_B0
        k_num = 2 * np.pi / L_num

        # RandomState(seed) draws the same stream so every rank gets identical phases 
        # global legacy rng is left untouched (pgen0 relies on that global rng to differ from rank to rank)
        phase_rng = np.random.RandomState(modes_perp)
        shape = (modes_perp, modes_perp, modes_par)
        ph1 = 2.0 * np.pi * phase_rng.rand(*shape)
        ph2 = 2.0 * np.pi * phase_rng.rand(*shape)
        ph3 = 2.0 * np.pi * phase_rng.rand(*shape)

        def perp_modes():
            # (bx, by) below is the curl of a scalar potential ~ cos(m kx) cos(n ky) sin(o kz),
            # so the perturbation is divergence free mode by mode.
            I = range(1, modes_perp + 1)
            for n, m in itertools.product(I, I):
                norm = np.sqrt(8.0) / (np.sqrt(n**2 + m**2) * modes_perp * np.sqrt(modes_par))
                for o in range(1, modes_par + 1):
                    yield (n, m, o, norm,
                           ph1[n - 1, m - 1, o - 1],
                           ph2[n - 1, m - 1, o - 1],
                           ph3[n - 1, m - 1, o - 1])

        def bx(x, y, z):
            b = np.zeros_like(x)
            for n, m, o, norm, p1, p2, p3 in perp_modes():
                b += norm * n * (np.sin(m * k_num * x + p1)
                                 * np.cos(n * k_num * y + p2)
                                 * np.sin(o * k_num * z + p3))
            return A0_num * b

        def by(x, y, z):
            b = np.zeros_like(x)
            for n, m, o, norm, p1, p2, p3 in perp_modes():
                b -= norm * m * (np.cos(m * k_num * x + p1)
                                 * np.sin(n * k_num * y + p2)
                                 * np.sin(o * k_num * z + p3))
            return A0_num * b

    # --------------------------------------------------
    # Print setup summary
    if runko.on_main_rank():
        from runko.auto_outdir import resolve_outdir
        W = 21
        logger.info(f"{'--- [io] ---':}")
        logger.info(f"  {'outdir':<{W}}= {conf.io_outdir}")
        logger.info(f"  {'resolved outdir':<{W}}= {resolve_outdir(conf)}")
        logger.info(f"  {'prefix / postfix':<{W}}= {conf.io_outdir_prefix} / {conf.io_outdir_postfix}")
        logger.info(f"  {'output_interval':<{W}}= {conf.io_output_interval}")
        logger.info(f"  {'spectra_nbins':<{W}}= {conf.io_n_spectra_bins}")
        logger.info(f"  {'spectra_umin':<{W}}= {conf.io_spectra_umin}")
        logger.info(f"  {'spectra_umax':<{W}}= {conf.io_spectra_umax}")

        logger.info(f"{'--- [grid] ---':}")
        logger.info(f"  {'tiles':<{W}}= {conf.n_tiles}")
        logger.info(f"  {'mesh per tile':<{W}}= {conf.n_cells_per_tile}")
        full_grid = np.array(conf.n_tiles) * np.array(conf.n_cells_per_tile)
        logger.info(f"  {'full grid':<{W}}= {full_grid}")

        logger.info(f"{'--- [simulation] ---':}")
        logger.info(f"  {'laps':<{W}}= {conf.n_laps}")
        logger.info(f"  {'eddy turnovers':<{W}}= {conf.n_eddy_turnovers}")
        logger.info(f"  {'eddy length':<{W}}= {eddy_length_num:.6g}")
        logger.info(f"  {'laps per turnover':<{W}}= {laps_per_eddy:.6g}")
        logger.info(f"  {'cfl':<{W}}= {conf.cfl}")

        logger.info(f"{'--- [particles] ---':}")
        logger.info(f"  {'ppc':<{W}}= {n_e_num}")
        logger.info(f"  {'m0 / m1':<{W}}= {conf.m0} / {conf.m1}")
        logger.info(f"  {'q0 / q1':<{W}}= {conf.q0:.6g} / {conf.q1:.6g}")
        logger.info(f"  {'theta_e / theta_i':<{W}}= {delgam0:.6g} / {delgam1:.6g}")

        logger.info(f"{'--- [problem] ---':}")
        logger.info(f"  {'sigma_e':<{W}}= {sigma_e}")
        logger.info(f"  {'sigma_i':<{W}}= {sigma_i}")
        logger.info(f"  {'sigma_e^hot':<{W}}= {sigma_e_hot}")
        logger.info(f"  {'sigma_i^hot':<{W}}= {sigma_i_hot}")
        logger.info(f"  {'sigma':<{W}}= {sigma}")
        logger.info(f"  {'h_e':<{W}}= {enth_e}")
        logger.info(f"  {'h_i':<{W}}= {enth_i}")
        logger.info(f"  {'alfven velocity':<{W}}= {valf:.6g}")
        logger.info(f"  {'B_init':<{W}}= {B0_num:.6g}")
        logger.info(f"  {'n_filter_passes':<{W}}= {conf.n_filter_passes}")

        if conf.use_antenna:
            logger.info(f"  {'turbulence':<{W}}= langevin antenna")
            logger.info(f"  {'modes_perp / _par':<{W}}= {modes_perp} / {modes_par}")
            logger.info(f"  {'A':<{W}}= {A_num:.6g}")
            logger.info(f"  {'linear freq':<{W}}= {linear_freq_num:.6g}")
            logger.info(f"  {'decorr freq':<{W}}= {gamma0_num:.6g}")
        else: # decaying
            logger.info(f"  {'turbulence':<{W}}= decaying")
            logger.info(f"  {'modes_perp / _par':<{W}}= {modes_perp} / {modes_par}")
            logger.info(f"  {'dB_perp / B0':<{W}}= {conf.Bperp_to_B0}")
            logger.info(f"  {'A0':<{W}}= {A0_num:.6g}")

        logger.info(f"{'--- [algorithms] ---':}")
        logger.info(f"  {'field_propagator':<{W}}= {conf.field_propagator}")
        logger.info(f"  {'particle_pusher':<{W}}= {conf.particle_pusher}")
        logger.info(f"  {'field_interpolator':<{W}}= {conf.field_interpolator}")
        logger.info(f"  {'current_depositer':<{W}}= {conf.current_depositer}")
        logger.info(f"  {'current_filter':<{W}}= {conf.current_filter}")
        logger.info(f"  {'tile_partitioning':<{W}}= {conf.tile_partitioning}")


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
                               bx, by, bz,
                               zero_field, zero_field, zero_field)

            if conf.use_antenna:
                for mode in antenna_modes:
                    logger.info("registering antenna...")
                    tile.register_antenna(mode)
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
        if conf.use_antenna:
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
