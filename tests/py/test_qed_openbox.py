# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analytic tests of the open box of the qed tile: sources, sinks and adaptivity.

  6 pair balance     thermostatted leptons, Wien photon injection, slab escape: creation = annihilation
  7 adaptivity       weight invariance; the resolution function N_LP/W ~ x^a
 12 corona           heavy-light accumulation; truncated emission above a floor; the hard-state cascade"""

import unittest
import numpy as np
import runko
from test_qed_twobody import (WL, H, NC, rows, set_pop, iso_dirs, spread, mj_pop, wien_pop, mono_pop,
                              init_pairs, energies, thermostat, reset_photons, totals, mj_ann_average, rho_gg)
from test_qed_onebody import box_N1 as box_N1_field, state, empty

ph = runko.qed.physics
LAMC_CM = 3.8615926796e-11


def make(seed, dt, H_lamC=H, w_lep=WL, ell_B=0.0, x_min=0.0, channels=None, nbin=8, nsort=1, slack=0.1,
         a_res=0.0, x_ref=1e-2, kappa=0.1, acc_syn=3.0, smin=4, nsub_max=8):
    c = runko.Configuration(None)
    c.n_tiles = [1, 1, 1]
    c.n_cells_per_tile = [NC, NC, NC]
    c.cfl = dt * NC
    c.field_propagator = "fdtd2"
    c.q0, c.m0, c.w0 = -1, 1, w_lep
    c.q1, c.m1, c.w1 = 1, 1, w_lep
    c.q2, c.m2 = 0, 0
    c.particle_pusher = "boris"
    c.field_interpolator = "linear_1st"
    c.current_depositer = "zigzag_1st"
    c.qed_H_lamC = H_lamC
    c.qed_ell_B = ell_B
    c.qed_x_min = x_min
    c.qed_nbin, c.qed_nsort, c.qed_slack, c.qed_smin = nbin, nsort, slack, smin
    c.qed_a_res, c.qed_x_ref = a_res, x_ref
    c.qed_kappa_acc, c.qed_acc_syn, c.qed_nsub_max = kappa, acc_syn, nsub_max
    c.qed_seed = seed
    if channels:
        c.qed_channels = list(channels)
    t = runko.qed.threeD.Tile((0, 0, 0), c)
    zero = lambda x, y, z: (0, 0, 0)
    t.set_EBJ(zero, zero, zero)
    return t


class qed_openbox(unittest.TestCase):

    def check(self, name, pred, meas, tol):
        err = abs(meas - pred) / (abs(pred) if pred != 0 else 1.0)
        self.assertLessEqual(err, tol, f"{name}: predicted {pred:.5g}, measured {meas:.5g}, tol {tol}")

    # ================================================================== test 6: pair balance
    def test6_pair_balance(self):
        """thermostatted leptons at Theta, Wien(Theta) photons injected at the compactness ell, slab escape;
        steady state: creation (spectrum integral) = annihilation (MJ rate) fixes W_pm and tau_T"""
        # the pair depth relaxes from 1.5 to its balance in ~2000 steps at the lowest compactness
        Theta, wl, Ne0, dt, steps = 0.5, 1.5e-4, 5000, 1 / 300, 5000
        N2 = 1 / dt
        bins = np.logspace(-2.5, 1.5, 45); cx = np.sqrt(bins[1:] * bins[:-1])
        rho_tab = np.array([[rho_gg(a, b) if a * b > 1 else 0.0 for b in cx] for a in cx])
        rho_ann = mj_ann_average(Theta, np.random.default_rng(7), 400000)
        for il, ell in enumerate((1.35, 4.05, 12.15)):
            rng = np.random.default_rng(52 + il)
            t = make(1, dt, w_lep=wl)
            um, _ = mj_pop(Ne0, Theta, rng); up, _ = mj_pop(Ne0, Theta, rng)
            set_pop(t, "e-", um, wl * np.ones(len(um))); set_pop(t, "e+", up, wl * np.ones(len(up)))
            ni = ell * dt / (3 * Theta * wl)          # injected photons per step at the weight w_lep
            uph, _ = wien_pop(int(ni * 300), Theta, rng); set_pop(t, "ph", uph, wl * np.ones(len(uph)))
            tau_rec, cre, ann, wpm, snap = [], [], [], [], []
            for n in range(steps + 1):
                tt = {s: t.totals(s) for s in ("e-", "e+")}
                tau_rec.append(tt["e-"]["W"] + tt["e+"]["W"]); wpm.append(tt["e-"]["W"] * tt["e+"]["W"])
                if n == steps:
                    break
                wm0 = tt["e-"]["W"]; c_ann = t.counts()["pairann"]
                thermostat(t, "e-", Theta, rng); thermostat(t, "e+", Theta, rng)
                t.interact()
                ann.append((t.counts()["pairann"] - c_ann) * wl)
                cre.append(t.totals("e-")["W"] - wm0 + ann[-1])   # created pair weight = dW_- + annihilated
                t.escape_photons()
                t.inject_photons(ell, Theta, wl, planck=False)
                if n >= 0.75 * steps and n % 25 == 0:
                    x, w = energies(t, "ph")
                    snap.append(np.histogram(x, bins=bins, weights=w)[0])
            late = slice(int(0.75 * steps), None)
            tau_m, cre_m, ann_m = np.mean(tau_rec[late]), np.mean(cre[late]), np.mean(ann[late])
            W_bin = np.mean(snap, 0)
            cre_pred = 0.5 * np.einsum("a,ab,b->", W_bin, rho_tab, W_bin) / N2
            ann_pred = np.mean(wpm[late]) * rho_ann / N2
            tau_pred = 2 * np.sqrt(cre_pred * N2 / rho_ann)   # pair balance from the measured photon spectrum
            self.check(f"tau_T/balance ell={ell:.1f}", 1.0, tau_m / tau_pred, 0.10)
            self.check(f"creation/integral ell={ell:.1f}", 1.0, cre_m / max(cre_pred, 1e-30), 0.10)
            self.check(f"annihilation/MJ ell={ell:.1f}", 1.0, ann_m / max(ann_pred, 1e-30), 0.10)

    def test_escape_external_depth(self):
        """the driver may give the slab's depth and step: the same numbers as the tile's own give the
        same escape, and the direction histogram sums to the spectrum"""
        dt = 5e-3
        rng = np.random.default_rng(21)
        tiles = [make(11, dt), make(11, dt)]
        for t in tiles:
            init_pairs(t, 4000, 1000, 0.1, 1.0, np.random.default_rng(21))
            t.sort_particles()
        tau, N2 = tiles[0].tau_T(), 1 / dt
        tiles[0].escape_photons()
        tiles[1].escape_photons(tau, dt)
        b0, b1 = tiles[0].budget(), tiles[1].budget()
        self.assertGreater(b0["esc_W"], 0.0)
        self.assertEqual(b0["esc_W"], b1["esc_W"])
        self.assertEqual(b0["esc_E"], b1["esc_E"])
        edges, h = tiles[0].escaped_spectrum()
        hmu = tiles[0].escaped_mu_spectrum()
        self.assertEqual(hmu.shape, (len(h), 32))
        np.testing.assert_allclose(hmu.sum(1), h, rtol=1e-12)
        self.assertAlmostEqual(h.sum(), b0["esc_E"], delta=1e-9 * b0["esc_E"])
        # a smaller step in crossings escapes proportionally fewer
        t = make(11, dt); init_pairs(t, 4000, 1000, 0.1, 1.0, np.random.default_rng(21)); t.sort_particles()
        t.escape_photons(tau, dt / 4)
        self.assertLess(t.budget()["esc_W"], 0.5 * b0["esc_W"])
        t.reset_budget()
        self.assertEqual(t.budget()["esc_W"], 0.0)
        self.assertEqual(t.escaped_mu_spectrum().sum(), 0.0)

    # ================================================================== test 7: adaptivity
    def test7_adaptivity(self):
        steps, dt = 300, 5e-3
        curves = {}
        for k in (1, 3, 10):
            # k times the weight on k times fewer LPs: the same physical box, hence the same evolution
            rng = np.random.default_rng(1)
            t = make(1, dt, w_lep=k * WL)
            u, w = wien_pop(8000 // k, 0.05, rng); set_pop(t, "ph", u, k * w)
            for s in ("e-", "e+"):
                u, w = mj_pop(2000 // k, 2.0, rng); set_pop(t, s, u, k * w)
            tt, E0, _ = totals(t)
            Eg = []
            for n in range(steps + 1):
                Eg.append(t.totals("ph")["E"] / E0)
                if n < steps:
                    t.interact()
            curves[k] = np.array(Eg)
        late = slice(int(0.75 * steps), None)
        ref = curves[1][late].mean()
        for k in (3, 10):
            self.check(f"weight invariance E_gamma k={k}", 1.0, curves[k][late].mean() / ref, 0.05)
        # (b) resolution function: N_LP/W relaxes to x^a once the roulette has thinned the soft end
        # and the copies have populated the hard one
        a, xref, steps_b = 0.3, 1e-2, 1500
        bins = np.logspace(-2, 1.3, 25); cx = np.sqrt(bins[1:] * bins[:-1])
        rng = np.random.default_rng(1)
        t = make(1, dt, a_res=a, x_ref=xref)
        init_pairs(t, 8000, 2000, 0.05, 2.0, rng)
        hN = np.zeros(len(cx)); hW = np.zeros(len(cx))
        for n in range(steps_b + 1):
            if n >= 0.75 * steps_b and n % 10 == 0:
                x, w = energies(t, "ph")
                hN += np.histogram(x, bins=bins)[0]; hW += np.histogram(x, bins=bins, weights=w)[0]
            if n < steps_b:
                t.interact()
        ok = hW > 50 * WL; r = hN[ok] / hW[ok]; xs = cx[ok]
        fit = np.polyfit(np.log(xs), np.log(r), 1)
        self.check("resolution slope a", a, fit[0], 0.15)

    # ================================================================== test 12: corona ingredients and the hard state
    def test12a_accumulation(self):
        """soft photons a thousand times heavier than the leptons: merged (kappa = 0.1) and plain (kappa = 0)
        events must move the same photon weight and the same lepton energy in expectation"""
        dt, Nph, Ne, x0, th, nseed, nstep = 1e-4, 2000, 1000, 1e-5, 1.0, 8, 100
        out = {}
        for kappa in (0.1, 0.0):
            wsc, dEl, nev = [], [], []
            for sd in range(nseed):
                rng = np.random.default_rng(62 + sd)
                t = make(1 + 100 * sd, dt, channels=["compton-"], kappa=kappa)
                set_pop(t, "ph", x0 * iso_dirs(Nph, rng), 1e3 * WL * np.ones(Nph))
                ue = np.array([0.0]); 
                from test_qed_twobody import sample_mj
                p = sample_mj(th, Ne, rng); set_pop(t, "e-", p[:, None] * iso_dirs(Ne, rng), WL * np.ones(Ne)); empty(t, ["e+"])
                dE = 0.0
                for n in range(nstep):
                    reset_photons(t, x0, rng); thermostat(t, "e-", th, rng)
                    E0 = t.totals("e-")["E"]; t.interact(); dE += t.totals("e-")["E"] - E0
                wsc.append(t.moved()["compton-"] / nstep); dEl.append(dE / nstep); nev.append(t.counts()["compton-"] / nstep)
            out[kappa] = (np.mean(wsc), np.std(wsc) / np.sqrt(nseed), np.mean(dEl), np.std(dEl) / np.sqrt(nseed), np.mean(nev))
        rng = np.random.default_rng(62)
        base = make(1, dt, channels=["compton-"], kappa=0.0)
        set_pop(base, "ph", x0 * iso_dirs(Nph, rng), 1e3 * WL * np.ones(Nph))
        from test_qed_twobody import sample_mj
        p = sample_mj(th, Ne, rng); set_pop(base, "e-", p[:, None] * iso_dirs(Ne, rng), WL * np.ones(Ne)); empty(base, ["e+"])
        exp_w = base.pair_sum("compton-", True)      # the weight the events move, independent of the merging
        m, p_ = out[0.1], out[0.0]
        self.assertLess(m[4], 0.5 * p_[4], "merging must cut the event count")
        self.check("accumulation: moved weight merged/plain", 1.0, m[0] / p_[0], max(0.03, 3 * (m[1] + p_[1]) / p_[0]))
        self.check("accumulation: moved weight/pair sum", 1.0, m[0] / exp_w, max(0.03, 3 * m[1] / exp_w))
        self.check("accumulation: lepton heating merged/plain", 1.0, m[2] / p_[2], max(0.05, 3 * (m[3] + p_[3]) / abs(p_[2])))

    def test12b_floor(self):
        """truncated emission: the emitted photons, their weight and their energy above a floor"""
        from test_qed_onebody import qed_onebody
        g0, N1, Ne2, nseed = 10.0, 0.02, 20000, 8
        b = float(ph.b_from_ellB(10.0, H)); chi = g0 * b
        xmin = 0.3 * 1.5 * g0 * g0 * b                # a floor below the critical energy: a third of the spectrum is cut
        rng = np.random.default_rng(72)
        pop = mono_pop(Ne2, g0, 1.0, rng)
        helper = qed_onebody()
        def box(seed):
            t = box_N1_field(seed, N1, ell_B=10.0, H=H, x_min=xmin, channels=["syn-"], acc_syn=0.0)
            set_pop(t, "e-", pop[0], pop[1]); empty(t, ["e+", "ph"])
            return t
        nex, wex, eex = helper.syn_expected(box(1), b, N1, x_min=xmin, acc=0.0)
        ev, wm, em, xs = [], [], [], []
        for sd in range(nseed):
            t = box(1 + 100 * sd); E0 = t.totals("e-")["E"]
            t.interact()
            ev.append(t.counts()["syn-"]); wm.append(t.totals("ph")["W"]); em.append(E0 - t.totals("e-")["E"])
            xs.append(state(t, "ph")[1])
        ev, wm, em = np.array(ev), np.array(wm), np.array(em)
        r0 = xmin / g0; P0 = float(ph.P_syn(chi, r0))
        x = np.sort(np.concatenate(xs)); r = np.clip(x / g0, r0, 1.0)
        cdf = (ph.P_syn(chi, r) - P0) / (1 - P0)
        ks = float(np.max(np.abs(cdf - np.arange(1, len(r) + 1) / len(r))))
        self.check("floor: events/expected", 1.0, ev.mean() / nex, 0.03)
        self.check("floor: emitted weight/expected", 1.0, wm.mean() / wex, 0.03)
        self.check("floor: energy lost/expected", 1.0, em.mean() / (Ne2 * WL) / eex, 0.03)
        self.check("floor: nothing emitted below it", 1.0, min(x.min() / xmin, 1.0), 1e-3)
        self.check("floor: spectrum KS vs truncated p_S", 0.0, ks, 4 / np.sqrt(len(r)))

    def test12c_corona(self):
        """the hard state of Nattila (2024): ell_th = 7, ell_nth = 1 (slope 6 on 3 < gamma < 30), ell_B = 0.8;
        a short run only checks that the box is heading to its steady state"""
        H_lamC, dt, n_cross = 9e6 / LAMC_CM, 1e-3, 3.0
        ell_th, ell_nth, p_inj, g_inj = 7.0, 1.0, 6.0, (3.0, 30.0)
        t = make(101, dt, H_lamC=H_lamC, ell_B=0.8, x_min=5e-6, nbin=16, nsort=10, slack=0.5, a_res=0.7, x_ref=1e-5)
        rng = np.random.default_rng(1)
        from test_qed_twobody import sample_mj
        Npm = int(round(0.5 * 1.0 / WL))
        for s in ("e-", "e+"):
            p = sample_mj(0.5, Npm, rng); set_pop(t, s, p[:, None] * iso_dirs(Npm, rng), WL * np.ones(Npm))
        empty(t, ["ph"])
        self.assertGreater(t.tau_ssa(5e-6), 1.0, "the floor must sit below the unit SSA depth")
        steps = int(round(n_cross / dt)); nlate = int(0.75 * steps)
        for n in range(steps):
            if n == nlate:
                t.reset_budget()
            t.interact()
            t.heat_diffusive(ell_th)
            t.inject_pairs_plaw(ell_nth, p_inj, g_inj[0], g_inj[1])
            t.escape_photons()
        bud = t.budget(); nl = steps - nlate
        ell_esc = bud["esc_E"] / nl / dt; ell_in = (bud["inj_E"] + bud["heat_E"]) / nl / dt
        self.check("corona: ell_esc/ell_inj", 1.0, ell_esc / (ell_th + ell_nth), 0.5)
        self.check("corona: injected power", 1.0, ell_in / (ell_th + ell_nth), 0.05)
        edges, h = t.escaped_spectrum(); cx = np.sqrt(edges[1:] * edges[:-1])
        ellx = h / nl / dt / np.diff(np.log(edges))
        x_peak = cx[np.argmax(ellx)]
        self.check("corona: log10 x_peak in [-0.5, 0.5]", 0.0, float(np.log10(x_peak)), 0.5)
        xg = np.geomspace(1e-7, 1e-2, 11); lt = np.log(np.maximum([t.tau_ssa(x) for x in xg], 1e-300))
        x_unit = np.exp(np.interp(0.0, lt[::-1], np.log(xg)[::-1]))
        self.check("corona: log10 x_SSA (unit SSA depth) in [-5.5, -4.5]", -5.0, float(np.log10(x_unit)), 0.1)
        self.assertGreater(t.tau_T(), 0.5)
        print(t.report())


if __name__ == "__main__":
    unittest.main()
