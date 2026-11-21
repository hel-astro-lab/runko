# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Radiative pair corona with the QED Monte Carlo of runko: the hard state of Nattila (2024).
Heated and injected pairs in a magnetized slab make every photon themselves; photons escape.
"""

import runko
import numpy as np
import argparse
import os
import time


if __name__ == "__main__":
    logger = runko.runko_logger()
    parser = argparse.ArgumentParser(description="QED Monte Carlo pair corona")
    parser.add_argument("--conf", type=str, default="corona.ini", help="Path to .ini configuration file")
    parser.add_argument("--n_cross", type=float, default=None, help="override the run length in crossings")
    args = parser.parse_args()

    conf = runko.Configuration(args.conf)
    if args.n_cross is not None:
        conf.n_cross = args.n_cross

    # --------------------------------------------------
    # The box: the tile edge in Compton wavelengths and the step in light crossings of it
    conf.qed_H_lamC = conf.H_cm / runko.qed.LAMC_CM
    dt = conf.cfl / conf.n_cells_per_tile[0]
    conf.n_laps = int(round(conf.n_cross / dt))
    w_lep = conf.w0
    ell_inj = conf.ell_th + conf.ell_nth
    rec_every = conf.io_record_interval or 10
    verbose_every = conf.io_verbose_interval or 500
    outdir = conf.io_outdir or "corona"

    if runko.on_main_rank():
        os.makedirs(outdir, exist_ok=True)
        W = 24
        logger.info("QED corona:")
        logger.info(f"  {'H':<{W}}= {conf.H_cm:.3e} cm = {conf.qed_H_lamC:.3e} lamC")
        logger.info(f"  {'dt':<{W}}= {dt:g} crossings ({conf.n_laps} laps for {conf.n_cross:g} crossings)")
        logger.info(f"  {'w_lep':<{W}}= {w_lep:g} Thomson depths")
        logger.info(f"  {'ell_th, ell_nth, ell_B':<{W}}= {conf.ell_th:g}, {conf.ell_nth:g}, {conf.qed_ell_B:g}")
        logger.info(f"  {'x_min':<{W}}= {conf.qed_x_min:g}")

    rng = np.random.default_rng(conf.qed_seed or 1)
    zero = lambda x, y, z: np.zeros_like(x)

    tile_grid = runko.TileGrid(conf)
    for idx in tile_grid.local_tile_indices():
        tile = runko.qed.threeD.Tile(idx, conf)
        tile.batch_set_EBJ(zero, zero, zero, zero, zero, zero, zero, zero, zero)
        # Maxwell-Juettner pairs at the initial depth tau0, no photons
        npm = int(round(0.5 * conf.tau0 / w_lep))
        for s in ("e-", "e+"):
            u = np.stack(runko.sample_boosted_juttner_synge(npm, conf.theta0, beta=0, gen=rng), 1)
            tile.set_particles(s, [tuple(v) for v in u], [w_lep] * npm)
        tile.set_particles("ph", [], [])
        b = tile.box()
        x_grid = np.geomspace(1e-7, 1e-2, 11)
        tss = [tile.tau_ssa(x) for x in [conf.qed_x_min] + list(x_grid)]
        logger.info(f"  N_2 = {b['N2']:.4g}, N_1 = {b['N1']:.4g}, b_perp = {b['b_perp']:.4e}")
        logger.info(f"  self-absorption depth of the initial leptons: tau_ssa(x_min = {conf.qed_x_min:.1e}) = {tss[0]:.3g}; "
                    + ", ".join(f"{a:.0e}: {c:.2g}" for a, c in zip(x_grid, tss[1:])))
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    # --------------------------------------------------
    # Records: the pair and energy budgets every rec_every laps; spectra over the last quarter
    rec = dict(t=[], tau=[], Wg=[], Eg=[], Epm=[], N=[], ell_esc=[], gbar=[], dtau_inj=[], counts=[], moved=[])
    prev = dict(counts=None, moved=None, esc_E=0.0, inj_W=0.0)
    lep_p, lep_w = [], []
    nlate = int(0.75 * conf.n_laps)
    t_start = time.time()

    def record(lap):
        for tile in simulation.local_tiles():
            # the state of the tile and its cumulative counters
            tt = {s: tile.totals(s) for s in ("ph", "e-", "e+")}
            bud = tile.budget(); cn = tile.counts(); mv = tile.moved()
            tau = tt["e-"]["W"] + tt["e+"]["W"]
            # the pair depth, the photon and lepton energies and weights, the LP counts
            rec["t"].append(lap * dt); rec["tau"].append(tau)
            rec["Wg"].append(tt["ph"]["W"]); rec["Eg"].append(tt["ph"]["E"]); rec["Epm"].append(tt["e-"]["E"] + tt["e+"]["E"])
            rec["gbar"].append((tt["e-"]["E"] + tt["e+"]["E"]) / max(tau, 1e-300))
            rec["N"].append([tt[s]["N"] for s in ("ph", "e-", "e+")])
            # rates since the previous record: escaping power, injected depth
            span = rec_every * dt
            rec["ell_esc"].append((bud["esc_E"] - prev["esc_E"]) / span if lap else 0.0)
            rec["dtau_inj"].append((bud["inj_W"] - prev["inj_W"]) / span if lap else 0.0)
            # events and moved weight per channel since the previous record
            rec["counts"].append({k: cn[k] - (prev["counts"][k] if prev["counts"] else 0.0) for k in cn})
            rec["moved"].append({k: mv[k] - (prev["moved"][k] if prev["moved"] else 0.0) for k in mv})
            prev.update(counts=cn, moved=mv, esc_E=bud["esc_E"], inj_W=bud["inj_W"])
            # lepton momenta over the last quarter, for the late distribution
            if lap >= nlate:
                vx, vy, vz = (np.asarray(a) for a in tile.get_velocities(tile.species("e-")))
                lep_p.append(np.sqrt(vx * vx + vy * vy + vz * vz)); lep_w.append(np.asarray(tile.get_weights(tile.species("e-"))))
                vx, vy, vz = (np.asarray(a) for a in tile.get_velocities(tile.species("e+")))
                lep_p.append(np.sqrt(vx * vx + vy * vy + vz * vz)); lep_w.append(np.asarray(tile.get_weights(tile.species("e+"))))

    def corona_step(x):
        lap = simulation.lap
        if lap % rec_every == 0:
            record(lap)
        if lap == nlate:
            for tile in simulation.local_tiles():
                tile.reset_budget()          # the escaping spectrum of the last quarter
            prev.update(esc_E=0.0, inj_W=0.0)
        # the reactions, then the sources and the sink of the open box
        x.prtcl_interact()
        x.grid_heat_diffusive(conf.ell_th)
        x.grid_inject_pairs_plaw(conf.ell_nth, conf.p_inj, conf.g_inj[0], conf.g_inj[1])
        x.grid_escape_photons()
        # the photons move; leptons too, in the fields (none here)
        x.prtcl_push()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)
        if lap % verbose_every == 0 and runko.on_main_rank():
            tile = next(iter(simulation.local_tiles()))
            tt = {s: tile.totals(s) for s in ("ph", "e-", "e+")}
            logger.info(f"lap {lap}/{conf.n_laps} t={lap*dt:.2f} tau_T={tt['e-']['W']+tt['e+']['W']:.3f} W_g={tt['ph']['W']:.3g} "
                        f"N_LP={[int(tt[s]['N']) for s in ('ph', 'e-', 'e+')]} ell_esc={rec['ell_esc'][-1]:.2f} ({time.time()-t_start:.0f}s)")

    simulation.for_each_lap(corona_step)
    record(conf.n_laps)

    # --------------------------------------------------
    # The steady state: energy balance, pair balance, spectral landmarks; saved for the plots
    tile = next(iter(simulation.local_tiles()))
    nl = conf.n_laps - nlate                 # laps of the last quarter
    # the escaping power and spectrum x ell_x per ln x, averaged over the last quarter
    bud = tile.budget()
    ell_esc = bud["esc_E"] / nl / dt
    edges, h = tile.escaped_spectrum(); cx = np.sqrt(edges[1:] * edges[:-1])
    ellx = h / nl / dt / np.diff(np.log(edges))
    # records as arrays; channel rates and the pair budget per crossing
    for k in ("t", "tau", "Wg", "Eg", "Epm", "N", "ell_esc", "gbar", "dtau_inj"):
        rec[k] = np.array(rec[k])
    names = list(rec["counts"][0])
    rates = {k: np.array([c[k] for c in rec["counts"]]) / (rec_every * dt) for k in names}
    dtau_cre = 2 * np.array([m.get("photann", 0.0) for m in rec["moved"]]) / (rec_every * dt)
    dtau_ann = 2 * np.array([m.get("pairann", 0.0) for m in rec["moved"]]) / (rec_every * dt)
    late = rec["t"] >= rec["t"][-1] * 0.75
    tau = float(rec["tau"][late].mean())
    # spectral landmarks: the peak, the soft turnover, the slope of the Comptonized band
    x_peak = cx[np.argmax(ellx)]
    soft = cx < 1e-3
    x_ssa = cx[soft][np.argmax(ellx[soft])] if soft.any() and ellx[soft].max() > 0 else np.nan
    band = (cx > 1e-3) & (cx < 1e-1) & (ellx > 0)
    slope = np.polyfit(np.log(cx[band]), np.log(ellx[band]), 1)[0] if band.sum() > 3 else np.nan
    # the late lepton distribution and its mean energy
    lp = np.concatenate(lep_p); lw = np.concatenate(lep_w)
    gbar = float(np.average(np.sqrt(1 + lp * lp), weights=lw))
    # the energy of unit self-absorption depth against the final leptons
    x_grid = np.geomspace(1e-7, 1e-2, 11)
    tss = np.array([tile.tau_ssa(x) for x in [conf.qed_x_min] + list(x_grid)])
    lt = np.log(np.maximum(tss[1:], 1e-300))
    x_unit = np.exp(np.interp(0.0, lt[::-1], np.log(x_grid)[::-1])) if lt.max() > 0 > lt.min() else np.nan
    # the photons inside the box at the end
    vx, vy, vz = (np.asarray(a) for a in tile.get_velocities(tile.species("ph")))
    ph_x = np.sqrt(vx * vx + vy * vy + vz * vz); ph_w = np.asarray(tile.get_weights(tile.species("ph")))
    inj, cre, ann = float(rec["dtau_inj"][late].mean()), float(dtau_cre[late].mean()), float(dtau_ann[late].mean())
    # the summary, the kernel report and the records for the plots
    if runko.on_main_rank():
        logger.info(f"{conf.n_laps} laps = {conf.n_cross:g} crossings in {time.time()-t_start:.0f}s; late average over the last quarter:")
        logger.info(f"   tau_T = {tau:.3f}; ell_esc = {ell_esc:.3f} against ell_inj = {ell_inj:g} (ratio {ell_esc/ell_inj:.3f}); "
                    f"W_g/tau_T = {rec['Wg'][late].mean()/tau:.3g}; eta_x = {rec['Eg'][late].mean()/tau:.3g}")
        logger.info(f"   pair balance, d tau_T/dt per crossing: injected {inj:.3f}, created by gamma-gamma {cre:.3f}, annihilated {ann:.3f} "
                    f"(sources/sink {(inj+cre)/max(ann, 1e-300):.3f})")
        logger.info(f"   escaping spectrum: peak of x ell_x at x = {x_peak:.3g}; soft turnover at x = {x_ssa:.3g}; "
                    f"slope over 1e-3 < x < 0.1: {slope:.2f} (photon index {2-slope:.2f})")
        logger.info(f"   self-absorption depth of the final leptons: tau_ssa(x_min) = {tss[0]:.3g}, unit depth at x = {x_unit:.2e}")
        logger.info(f"   leptons: <gamma> = {gbar:.3f}; events per crossing " + ", ".join(f"{k}: {rates[k][late].mean():.3g}" for k in names))
        logger.info(tile.report())
        np.savez_compressed(os.path.join(outdir, "corona.npz"),
                            t=rec["t"], tau=rec["tau"], Wg=rec["Wg"], Eg=rec["Eg"], Epm=rec["Epm"], N=rec["N"],
                            ell_esc=rec["ell_esc"], gbar=rec["gbar"], dtau_cre=dtau_cre, dtau_ann=dtau_ann, dtau_inj=rec["dtau_inj"],
                            rate_names=np.array(names), rates=np.array([rates[k] for k in names]),
                            cx=cx, ellx=ellx, lep_p=lp, lep_w=lw, ph_x=ph_x, ph_w=ph_w,
                            scalars=np.array([tau, ell_esc, ell_inj, x_peak, x_ssa, slope, gbar, nl, len(lep_p) // 2, w_lep, conf.qed_x_min, x_unit, tss[0], cre, ann, inj, conf.qed_a_res or 0.7, conf.qed_x_ref or 1e-5]))
        logger.info(f"saved {os.path.join(outdir, 'corona.npz')}")
