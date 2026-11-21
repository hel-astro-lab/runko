# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The pair corona inside the PIC loop: a periodic box driven by the Langevin antenna (as
projects/pic-turbulence/pic-langevin-antenna.py), the QED Monte Carlo of every channel acting on
the local field at the top of each lap, photons escaping through the slab sink. Energy enters
only through the antenna; the turbulence it drives energizes the pairs, which radiate.
"""

import runko
import numpy as np
import argparse
import os
import time
from mpi4py import MPI

from qed_diagnostics import QedDiagnostics


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    logger = runko.runko_logger()

    parser = argparse.ArgumentParser(description="QED Monte Carlo pair corona in PIC turbulence")
    parser.add_argument("--conf", type=str, default="corona-pic.ini", help="Path to .ini configuration file")
    parser.add_argument("--n_cross", type=float, default=None, help="override the run length in crossings")
    args = parser.parse_args()

    conf = runko.Configuration(args.conf)
    if args.n_cross is not None:
        conf.n_cross = args.n_cross
    # every rank draws the same antenna: a shared seed
    rng = np.random.default_rng(seed=conf.qed_seed or 1)

    # --------------------------------------------------
    # The box: the tile edge sets the QED clocks, the box edge the escape and the compactness
    n_t = conf.n_tiles[0]
    n_cells = conf.n_cells_per_tile[0]
    if conf.n_tiles != [n_t] * 3 or conf.n_cells_per_tile != [n_cells] * 3:
        raise RuntimeError("qed-pic needs a cubic grid of cubic tiles")
    conf.qed_H_lamC = conf.H_cm / n_t / runko.qed.LAMC_CM   # tile edge in Compton wavelengths
    dt_tile = conf.cfl / n_cells                     # step in crossings of a tile
    dt_box = dt_tile / n_t                           # and of the box
    conf.n_laps = int(round(conf.n_cross / dt_box))
    # lepton weight from the initial depth: tau_box = sum w / n_t^2 over both species
    n_cells_total = (n_t * n_cells) ** 3
    n_lep = 2 * conf.ppc * n_cells_total
    w_lep = conf.tau0 * n_t ** 2 / n_lep
    conf.w0 = conf.w1 = w_lep
    ell_B = conf.qed_ell_B

    # --------------------------------------------------
    # The plasma and the antenna, as in the PIC turbulence driver
    n_e_num = conf.ppc
    modes_par, modes_perp = conf.modes_par, conf.modes_perp
    modes = [(modes_perp, 0, modes_par), (modes_perp, 0, -modes_par), (0, modes_perp, modes_par), (0, modes_perp, -modes_par)]
    delgam = conf.theta0
    temp_ratio = conf.theta1_to_theta0 or 1.0
    n_num = 2 * n_e_num
    plasma_freq_num = conf.cfl / conf.n_cells_per_skindepth
    average_gamma = 1 + 3 * delgam
    conf.q0 = -average_gamma * (plasma_freq_num ** 2.0) / (n_e_num * (1.0 + conf.m0 / conf.m1))
    conf.q1 = abs(conf.q0)
    m0_num = conf.m0 * abs(conf.q0)
    delgam0 = delgam
    delgam1 = temp_ratio * delgam0
    sigma = conf.sigma / average_gamma
    B0_num = np.sqrt(n_num * (conf.cfl ** 2.0) * sigma * m0_num)
    conf.qed_B0 = B0_num                             # the field of ell_B, in code units
    L_num = n_t * n_cells
    kpar_num = 2 * np.pi * modes_par / L_num
    kperp_num = 2 * np.pi * modes_perp / L_num
    A_num = kpar_num * B0_num / (kperp_num ** 2 * np.sqrt(len(modes)))
    valf = np.sqrt(sigma / (1 + sigma))
    linear_freq_num = kpar_num * conf.cfl * valf
    w0_num = 0.8 * linear_freq_num
    gamma0_num = -0.6 * linear_freq_num
    base_time_evolution = runko.sample_oscillating_langevin_antenna(size=conf.n_laps, characteristic_freq=w0_num,
                                                                    decorrelation_rate=gamma0_num, gen=rng)

    def antenna(mode):
        random_phase = np.exp(2j * np.pi * rng.random())
        if conf.fixed_antenna:
            return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=base_time_evolution * random_phase)
        time_evolution = runko.sample_oscillating_langevin_antenna(size=conf.n_laps, characteristic_freq=w0_num,
                                                                   decorrelation_rate=gamma0_num, gen=rng)
        return runko.emf.threeD.antenna_mode(A=(0, 0, A_num), n=mode, lap_coeffs=time_evolution * random_phase)

    antenna_modes = list(map(antenna, modes))

    n_filter_passes = conf.n_filter_passes or 3
    output_interval = conf.io_output_interval or 100
    qed_interval = conf.io_qed_interval or 20
    verbose_interval = conf.io_verbose_interval or 100
    from runko.auto_outdir import resolve_outdir
    outdir = resolve_outdir(conf)

    if runko.on_main_rank():
        os.makedirs(outdir, exist_ok=True)
        W = 24
        logger.info("QED corona in PIC turbulence:")
        logger.info(f"  {'box':<{W}}= {n_t}^3 tiles of {n_cells}^3 cells, H = {conf.H_cm:.3e} cm")
        logger.info(f"  {'dt':<{W}}= {dt_box:.3e} box crossings ({conf.n_laps} laps for {conf.n_cross:g} crossings)")
        logger.info(f"  {'w_lep':<{W}}= {w_lep:.3e} (tau0 = {conf.tau0:g} over {n_lep} leptons)")
        logger.info(f"  {'ell_B, x_min':<{W}}= {ell_B:g}, {conf.qed_x_min:g}")
        logger.info(f"  {'sigma, B0, A':<{W}}= {conf.sigma:g}, {B0_num:.4g}, {A_num:.4g}")
        logger.info(f"  {'q0, m0':<{W}}= {conf.q0:.4g}, {conf.m0}")
        logger.info(f"  {'outdir':<{W}}= {outdir}")

    # --------------------------------------------------
    # Fields and particles: the mean field along z, co-located thermal pairs, no photons
    zero_field = lambda x, y, z: np.zeros_like(x)
    bz = lambda x, y, z: np.ones_like(x) * B0_num

    def pgen0(x, y, z):
        N = len(x)
        # different ranks draw different positions: bare np.random, not the shared rng
        pgen0.pos = x + np.random.random(N), y + np.random.random(N), z + np.random.random(N)
        vel = runko.sample_boosted_juttner_synge(N, delgam0, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    def pgen1(x, y, z):
        vel = runko.sample_boosted_juttner_synge(len(x), delgam1, beta=0)
        return runko.pic.threeD.ParticleStateBatch(pos=pgen0.pos, vel=vel)

    tile_grid = runko.TileGrid(conf)
    for idx in tile_grid.local_tile_indices():
        tile = runko.qed.threeD.Tile(idx, conf)
        tile.batch_set_EBJ(zero_field, zero_field, zero_field, zero_field, zero_field, bz, zero_field, zero_field, zero_field)
        for a in antenna_modes:
            tile.register_antenna(a)
        for _ in range(n_e_num):
            tile.batch_inject_to_cells(0, pgen0)
            tile.batch_inject_to_cells(1, pgen1)
        tile.set_particles("ph", [], [])
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    def sync_EB(x):
        EB = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)
        x.comm_external(*EB)
        x.comm_local(*EB)

    simulation.prelude(sync_EB)

    diag = QedDiagnostics(simulation, dt_box, n_t, outdir)
    t_start = time.time()
    last_record = 0

    def global_tau():
        """the pair depth of the whole box: the lepton weight summed over every rank"""
        loc = sum(tile.totals("e-")["W"] + tile.totals("e+")["W"] for tile in simulation.local_tiles())
        return comm.allreduce(loc, op=MPI.SUM) / n_t ** 2

    def corona_step(x):
        lap = simulation.lap
        # ---- QED first: every particle is inside its tile after the last exchange, products land
        # at their parents' positions and take the push and the exchange of this lap
        x.prtcl_interact()
        x.grid_escape_photons(global_tau(), dt_box)
        # ---- the PIC lap of the antenna-driven turbulence
        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        x.prtcl_push()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

        x.prtcl_deposit_current()
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J_exchange)
        x.comm_external(runko.tools.comm_mode.emf_J)
        x.comm_local(runko.tools.comm_mode.emf_J)
        for i in range(n_filter_passes):
            if i > 0 and i % 3 == 0:
                x.comm_external(runko.tools.comm_mode.emf_J)
                x.comm_local(runko.tools.comm_mode.emf_J)
            x.grid_filter_current()

        x.grid_push_half_b()
        x.comm_external(runko.tools.comm_mode.emf_B)
        x.comm_local(runko.tools.comm_mode.emf_B)

        x.grid_push_e()
        x.grid_deposit_antenna_current()
        x.grid_add_current()
        x.comm_external(runko.tools.comm_mode.emf_E)
        x.comm_local(runko.tools.comm_mode.emf_E)

        # ---- diagnostics and output
        x.io_average_kinetic_energy()
        x.io_average_B_energy_density()
        x.io_average_E_energy_density()
        global last_record
        if lap % qed_interval == 0:
            diag.record(lap, lap - last_record)
            last_record = lap
        if lap % output_interval == 0:
            x.io_emf_snapshot()
            x.io_prtcl_snapshot()
            x.io_spectra_snapshot()
        if lap % verbose_interval == 0 and runko.on_main_rank() and diag.series["tau"]:
            S = diag.series
            logger.info(f"lap {lap}/{conf.n_laps} t={lap*dt_box:.3f} tau_T={S['tau'][-1]:.3f} [{S['tau_min'][-1]:.2f}, {S['tau_max'][-1]:.2f}] "
                        f"N_LP={[int(n) for n in S['N'][-1]]} ell_esc={S['ell_esc'][-1]:.3g} ({time.time()-t_start:.0f}s)")

    simulation.for_each_lap(corona_step)
    # the laps 0 .. n_laps - 1 ran: n_laps - 1 - last_record of them since the last record
    if conf.n_laps - 1 > last_record:
        diag.record(conf.n_laps, conf.n_laps - 1 - last_record)
    diag.save_series()
    if runko.on_main_rank():
        tile = next(iter(simulation.local_tiles()))
        logger.info(tile.report())
        logger.info(f"{conf.n_laps} laps in {time.time()-t_start:.0f}s; output in {outdir}")
