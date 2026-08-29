# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analytic tests of the field channels of the qed tile.

  9 self-absorption      the 2->1 channel: rate vs the classical form and the pair sum; 2->1 bookkeeping;
                         the bin bound covers every proposable pair
 10 synchrotron          the one-body emission channel: event rate, accumulation, cooling, energy budget
 11 magnetic pair prod.  the one-body conversion channel: event rate, 1->2 bookkeeping, charge symmetry
The field enters through the magnetic compactness and the tile size; the one-body clock is
N_1 = (lamC/H)/(alpha dt), so a box is built from N_1 by choosing its step."""

import unittest
import numpy as np
from scipy.special import kv
import runko

from tests.qed_test_helpers import (H, WL, make_tile, box_N1, set_pop, empty, iso_dirs, spread, mono_pop, state,
                                    totals, CheckMixin)


ph = runko.qed.physics


class qed_onebody(CheckMixin, unittest.TestCase):

    # ================================================================== test 9: self-absorption
    def test9_self_absorption(self):
        dt, ell_B, x_min = 2e-3, 10.0, 1e-2
        g0, Nph, Ne, nseed = 10.0, 4000, 1500, 16
        b = float(ph.b_from_ellB(ell_B, H))
        # ---- the rate against the classical closed form, at gammas where its 1/gamma^2 correction is
        # below the accuracy of the fitted synchrotron functions
        bc = 1e-10
        for g in (30.0, 100.0, 300.0):
            nu = np.geomspace(1e-2, 3.0, 60); x = 1.5 * g * g * bc * nu
            num = ph.rate_ssa(x, 0 * x, 0 * x, np.sqrt(g * g - 1), 0.0, 0.0, bc, 0.0)
            cls = (3 * np.sqrt(3) / (16 * ph.ALPHA)) * bc * nu ** 2 * kv(5 / 3, nu) / (x * x * g)
            self.check(f"rate vs classical, gamma={g:.0f}", 0.0, float(np.max(np.abs(num / cls - 1))), 0.01)
        # ---- absorbed weight against the exact pair sum: merged events move f s each
        rng = np.random.default_rng(92)
        ue, we = mono_pop(Ne, g0, 1.0, rng)
        xc = 1.5 * g0 * g0 * b
        x = xc * np.exp(rng.uniform(np.log(0.05), np.log(5.0), Nph))
        uph, wph = x[:, None] * iso_dirs(Nph, rng), spread(Nph, rng)
        base = make_tile(1, dt, ell_B=ell_B, x_min=x_min, channels=["ssa-"])
        set_pop(base, "ph", uph, wph); set_pop(base, "e-", ue, we); empty(base, ["e+"])
        exp = base.pair_sum("ssa-", True)
        meas, nev = [], []
        for sd in range(nseed):
            t = make_tile(1 + 100 * sd, dt, ell_B=ell_B, x_min=x_min, channels=["ssa-"])
            set_pop(t, "ph", uph, wph); set_pop(t, "e-", ue, we); empty(t, ["e+"])
            t.interact(); meas.append(t.moved()["ssa-"]); nev.append(t.counts()["ssa-"])
        meas = np.array(meas)
        self.check(f"absorbed weight/expected (first step, +- {meas.std()/np.sqrt(nseed)/exp:.3f})", 1.0, meas.mean() / exp, 0.04)
        self.assertGreater(meas.mean() / np.mean(nev) / WL, 1.0)  # merged: more weight than events x w_lep
        # ---- a long run: energy transfer, weight bookkeeping
        steps = 60
        t = make_tile(1, dt, ell_B=ell_B, x_min=x_min, channels=["ssa-"])
        set_pop(t, "ph", uph, wph); set_pop(t, "e-", ue, we); empty(t, ["e+"])
        t0, E0 = totals(t)
        for n in range(steps):
            t.interact()
        t1, E1 = totals(t)
        self.check("total energy conserved", 0.0, abs(E1 - E0) / E0, 1e-3)
        self.check("lepton weight conserved", 0.0, abs(t1["e-"]["W"] - t0["e-"]["W"]) / t0["e-"]["W"], 1e-12)
        loss = t0["ph"]["E"] - t1["ph"]["E"]
        self.assertGreater(loss, 0.05 * t0["ph"]["E"])
        self.check("photon energy loss = lepton gain", 0.0, abs(loss - (t1["e-"]["E"] - t0["e-"]["E"])) / loss, 0.05)
        # ---- the bound must cover the rate over every pair the envelope can propose, bin by bin
        t.sort_particles()
        tp, te = t.get_tables("ph"), t.get_tables("e-")
        kp, ke = np.asarray(t.get_slot_keys("ph")), np.asarray(t.get_slot_keys("e-"))
        up, xp, _ = state(t, "ph"); ue1, ge, _ = state(t, "e-")
        nb = 8
        rg = np.random.default_rng(8)
        ii = rg.integers(0, len(xp), 200000); jj = rg.integers(0, len(ge), 200000)
        r = ph.rate_ssa(up[ii, 0], up[ii, 1], up[ii, 2], ue1[jj, 0], ue1[jj, 1], ue1[jj, 2], b, x_min)
        a = (kp[ii] // 2) % nb; bb = ke[jj] % nb
        # the bound of the (photon bin a, lepton bin b) pair from the union of the two photon bits
        xlo = np.minimum(tp["emin"][2 * a], tp["emin"][2 * a + 1])
        sh = ph.ssa_bound(b, x_min, xlo, te["emin"][bb], te["emax"][bb])
        worst = float(np.max(r / np.maximum(sh, 1e-300)))
        self.check("max rate / bound <= 1", 0.0, max(worst - 1.0, 0.0), 1e-12)

    # ================================================================== test 10: quantum synchrotron
    def syn_expected(self, t, b, N1, x_min=0.0, acc=3.0, kappa=0.1):
        """expected events, emitted weight and energy loss of the emission channel in one step"""
        u, gam, w = state(t, "e-")
        lam = ph.rate_syn(u[:, 0], u[:, 1], u[:, 2], b, x_min, False) / N1
        xm = ph.xmean_syn(gam * b, gam, False, x_min)
        slow = lam <= 1.0
        f = np.where(slow, ph.facc_syn(gam, xm, acc, kappa), 1.0)
        p = np.where(slow, 1 - np.exp(-lam / f), 1.0)
        mult = np.where(slow, 1.0, (gam - ph.gamma_cooled(gam, lam, xm, b, N1)) / np.maximum(xm, 1e-300))
        return float(p.sum()), float(np.sum(w * f * p * mult)), float(np.mean(p * f * mult * xm))

    def syn_box(self, seed, N1, pop, acc=3.0):
        t = box_N1(seed, N1, channels=["syn-"], acc_syn=acc)
        set_pop(t, "e-", pop[0], pop[1]); empty(t, ["e+", "ph"])
        return t

    def test10_synchrotron(self):
        g0, Ne, nseed = 30.0, 20000, 8
        N1r = 0.02                                   # resolved: lambda ~ 0.07 per lepton and step
        b = float(ph.b_from_ellB(10.0, H)); chi = g0 * b
        rng = np.random.default_rng(32)
        pop = mono_pop(Ne, g0, 1.0, rng)
        lam0 = float(ph.rate_syn(np.sqrt(g0 * g0 - 1), 0.0, 0.0, b, 0.0, False)) / N1r
        t = self.syn_box(1, N1r, pop)
        self.assertAlmostEqual(t.box()["N1"], N1r, delta=1e-9 * N1r)
        # ---- rate x mean drawn photon = radiated power, times the drawn/exact moment ratio rho
        rho = float(ph.rmean_syn(chi, 0.0) * 3 * ph.T_syn(chi) / (2 * chi * ph.g_syn(chi)))
        pw = lam0 * float(ph.xmean_syn(chi, g0, False, 0.0)) / ((2 / 3) * chi * chi * float(ph.g_syn(chi)) / N1r)
        self.check("rate x mean photon = radiated power", rho, pw, 1e-9)
        # ---- event rate and emitted weight against the per-LP expectation, with and without merging
        nlp = {}
        for acc, tag in ((3.0, "merged"), (0.0, "plain")):
            nex, wex, _ = self.syn_expected(self.syn_box(1, N1r, pop, acc), b, N1r, acc=acc)
            ev, wm, nl = [], [], []
            for sd in range(nseed):
                t = self.syn_box(1 + 100 * sd, N1r, pop, acc); t.interact()
                ev.append(t.counts()["syn-"]); wm.append(t.totals("ph")["W"]); nl.append(t.totals("ph")["N"])
            ev, wm = np.array(ev), np.array(wm)
            self.check(f"events/expected ({tag})", 1.0, ev.mean() / nex, max(0.02, 4 * ev.std() / np.sqrt(nseed) / nex))
            self.check(f"emitted weight/expected ({tag})", 1.0, wm.mean() / wex, max(0.02, 4 * wm.std() / np.sqrt(nseed) / wex))
            nlp[tag] = np.mean(nl)
        xm0 = float(ph.xmean_syn(chi, g0, False, 0.0))
        f_pred = float(ph.facc_syn(g0, xm0, 3.0, 0.1))
        self.check("LP saving / f_acc", 1.0, nlp["plain"] / max(nlp["merged"], 1) / f_pred, 0.06)
        # ---- a long resolved run: the parent gives up exactly what the photons carry
        t = self.syn_box(1, N1r, pop)
        t0, E0 = totals(t)
        for _ in range(60):
            t.interact()
        t1, E1 = totals(t)
        self.check("total energy conserved (resolved)", 0.0, abs(E1 - E0) / E0, 1e-6)
        self.check("lepton loss = photon gain", 0.0, abs((t0["e-"]["E"] - t1["e-"]["E"]) - t1["ph"]["E"]) / t1["ph"]["E"], 1e-6)
        self.check("lepton weight conserved", 0.0, abs(t1["e-"]["W"] - t0["e-"]["W"]) / t0["e-"]["W"], 1e-12)
        # ---- the drawn spectrum: one step of a large population at a single chi_e
        t = self.syn_box(6, 0.003, pop, acc=0.0)      # lambda ~ 0.47, still resolved
        t.interact()
        _, x, _ = state(t, "ph")
        r = np.sort(np.clip(x / g0, 1e-12, 1.0))
        ks = float(np.max(np.abs(ph.P_syn(chi, r) - np.arange(1, len(r) + 1) / len(r))))
        self.check("drawn spectrum KS vs p_S", 0.0, ks, 4 / np.sqrt(len(r)))
        # ---- the unresolved regime: the parent must follow the closed cooling law
        N1u, steps_u = 1e-4, 30                       # lambda ~ 14 per lepton and step
        t = self.syn_box(8, N1u, pop)
        gh = [t.totals("e-")["E"] / t.totals("e-")["W"]]
        for _ in range(steps_u):
            t.interact(); gh.append(t.totals("e-")["E"] / t.totals("e-")["W"])
        gp = [g0]
        for _ in range(steps_u):
            g = gp[-1]
            if g >= ph.SYN_GMIN:
                lam = float(ph.rate_syn(np.sqrt(g * g - 1), 0.0, 0.0, b, 0.0, False)) / N1u
                xm = float(ph.xmean_syn(g * b, g, False, 0.0))
                g = float(ph.gamma_cooled(g, lam, xm, b, N1u))
            gp.append(g)
        gh, gp = np.array(gh), np.array(gp)
        self.check("unresolved gamma vs cooling law", 0.0, float(np.max(np.abs(gh / gp - 1))), 1e-5)
        t1, _ = totals(t)
        self.check("unresolved loss = gain (in expectation)", 0.0, abs((g0 * t1["e-"]["W"] - t1["e-"]["E"]) - t1["ph"]["E"]) / t1["ph"]["E"], 0.02)
        # ---- across the transition: the energy given up in one step against the two closed forms
        dev = 0.0
        for N1 in np.geomspace(0.015, 1e-5, 12):
            t = self.syn_box(14, N1, pop, acc=0.0)
            g_before = t.totals("e-")["E"] / t.totals("e-")["W"]
            _, _, pr = self.syn_expected(t, b, N1, acc=0.0)
            t.interact()
            dev = max(dev, abs((g_before - t.totals("e-")["E"] / t.totals("e-")["W"]) / pr - 1))
        self.check("loss/predicted across lambda = 1", 0.0, dev, 0.10)

    # ================================================================== test 11: magnetic pair production
    def bw_box(self, seed, N1, x0, Nph, spread_x=1.0):
        t = box_N1(seed, N1, ell_B=2.4323, H_lamC=1e5, channels=["bw"])
        rng = np.random.default_rng(seed + 17)
        x = x0 * (np.exp(rng.uniform(-np.log(spread_x), np.log(spread_x), Nph)) if spread_x > 1 else np.ones(Nph))
        set_pop(t, "ph", x[:, None] * iso_dirs(Nph, rng), spread(Nph, rng)); empty(t, ["e-", "e+"])
        return t

    def test11_magnetic_pairs(self):
        N1, x0, Nph, nseed = 0.5, 20.0, 20000, 8    # b_perp = 0.1: chi_gamma = 2, lambda ~ 1e-2 per photon and step
        b = float(ph.b_from_ellB(2.4323, 1e5)); chi = x0 * b
        self.assertAlmostEqual(b, 0.1, places=4)
        xb = np.geomspace(1e-2, 1.999, 50)
        self.check("rate is zero below threshold", 0.0, float(np.max(ph.rate_bw(xb, 0 * xb, 0 * xb, b))), 1e-300)
        base = self.bw_box(1, N1, x0, Nph, 3.0)
        u, x, _ = state(base, "ph")
        exp = float(np.sum(1 - np.exp(-ph.rate_bw(u[:, 0], u[:, 1], u[:, 2], b) / N1)))
        ev = []
        for sd in range(nseed):
            t = self.bw_box(1 + 100 * sd, N1, x0, Nph, 3.0); t.interact(); ev.append(t.counts()["bw"])
        ev = np.array(ev)
        self.check("conversions/expected (first step)", 1.0, ev.mean() / exp, max(0.05, 4 * ev.std() / np.sqrt(nseed) / exp))
        # ---- a long run: 1 -> 2 bookkeeping, and the photon decay against exp(-lambda n)
        steps = 60
        t = self.bw_box(1, N1, x0, Nph)
        lam0 = float(ph.rate_bw(x0, 0.0, 0.0, b)) / N1
        t0, E0 = totals(t); hist = []
        for n in range(steps + 1):
            hist.append(t.totals("ph")["W"] / t0["ph"]["W"])
            if n < steps:
                t.interact()
        hist = np.array(hist); t1, E1 = totals(t)
        lossE = t0["ph"]["E"] - t1["ph"]["E"]; lossW = t0["ph"]["W"] - t1["ph"]["W"]
        tol = 2 * WL / lossW    # the pair count is rounded systematically per cell: one lepton per cell
        self.check("total energy conserved", 0.0, abs(E1 - E0) / E0, tol)
        self.check("photon energy loss = pair gain", 0.0, abs(lossE - (t1["e-"]["E"] + t1["e+"]["E"])) / lossE, tol)
        self.check("photon weight loss = pair weight", 0.0, abs(lossW - t1["e-"]["W"]) / lossW, tol)
        self.check("charge symmetry W_- = W_+", 0.0, abs(t1["e-"]["W"] - t1["e+"]["W"]), 1e-12)
        self.check("photon decay vs exp(-lambda n)", 0.0, float(np.max(np.abs(hist / np.exp(-lam0 * np.arange(steps + 1)) - 1))), 0.05)
        # ---- the drawn pair spectrum: one step of a large population at a single chi_gamma
        t = self.bw_box(6, 0.01, x0, 20000)          # lambda ~ 0.5
        t.interact()
        _, gm, wm = state(t, "e-"); _, gp, _ = state(t, "e+")
        o = np.argsort(gm); r = np.clip((gm[o] - 1) / (x0 - 2), 0.0, 1.0); cw = np.cumsum(wm[o]) / wm.sum()
        ks = float(np.max(np.abs(ph.P_bw(chi, r) - cw)))
        self.check("drawn spectrum KS vs p_BW", 0.0, ks, 4 / np.sqrt(len(r)))
        self.check("mean chi_-/chi_gamma = 1/2", 0.5, float(np.average(r, weights=wm[o])), 8 * float(r.std()) / np.sqrt(len(r)))
        a, bb = np.sort(gm), np.sort(gp); ab = np.concatenate([a, bb])
        sym = float(np.max(np.abs(np.searchsorted(a, ab, "right") / len(a) - np.searchsorted(bb, ab, "right") / len(bb))))
        self.check("e- and e+ drawn alike (two-sample KS)", 0.0, sym, 4 * np.sqrt(2 / len(r)))


if __name__ == "__main__":
    unittest.main()
