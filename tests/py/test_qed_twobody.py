# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analytic tests of the two-body QED channels of the qed tile (the closed box).

  1 boosted box          equilibrium from the invariant mass of a closed box, all four channels
  2 Compton thermostat   Wien relaxation; Kompaneets log-normal; amplification per scattering A(Theta)
  3 IC cooling           gamma(t) of monoenergetic electrons in a soft photon bath vs the KN-corrected law
  4 pair annihilation    cold decay law; annihilation line; thermal rate vs the MJ average
  5 photon-photon        two monoenergetic photon gases vs the rate ODE
  8 rate exactness       first-step event counts vs the exact pair sum, versus the step
Every test builds its box from physical inputs (the step in light crossings, the lepton weight in
Thomson depths) and states its prediction through the clock N_2 = 1/dt the box derives."""

import unittest
import numpy as np
from scipy.special import kv
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import runko

from tests.qed_test_helpers import (WL, make_tile, rows, set_pop, iso_dirs, mj_momenta, spread, mono_pop,
                                    CheckMixin)


# ======================================================================================= helpers
def N2(dt):
    return 1.0 / dt


def rows(u):
    return [tuple(map(float, v)) for v in np.asarray(u).reshape(-1, 3)]


def set_pop(t, s, u, w):
    t.set_particles(s, rows(u), [float(x) for x in np.asarray(w)])


def iso_pairs(n, rng):
    v = iso_dirs(n // 2, rng)
    return np.concatenate([v, -v])


def mj_pop(n, theta, rng):
    """Maxwell-Juettner leptons in momentum-balanced pairs"""
    n -= n % 2
    u = mj_momenta(n // 2, theta, rng)
    return np.concatenate([u, -u]), WL * np.ones(n)


def wien_pop(n, theta, rng):
    n -= n % 2
    x = rng.gamma(3.0, theta, n // 2); x = np.concatenate([x, x])
    return x[:, None] * iso_pairs(n, rng), spread(n, rng)


def init_pairs(t, Nph, Npm, th_ph, th_pm, rng):
    u, w = wien_pop(Nph, th_ph, rng); set_pop(t, "ph", u, w)
    for s in ("e-", "e+"):
        u, w = mj_pop(Npm, th_pm, rng); set_pop(t, s, u, w)


def energies(t, s):
    m = 0.0 if s == "ph" else 1.0
    vx, vy, vz = (np.asarray(a) for a in t.get_velocities(t.species(s)))
    return np.sqrt(m * m + vx * vx + vy * vy + vz * vz), np.asarray(t.get_weights(t.species(s)))


def thermostat(t, s, theta, rng):
    n = int(t.totals(s)["N"])
    t.reset_momenta(s, rows(mj_momenta(n, theta, rng)))


def reset_photons(t, x, rng):
    n = int(t.totals("ph")["N"])
    xx = np.asarray(x) if np.ndim(x) else x * np.ones(n)
    t.reset_momenta("ph", rows(xx[:, None] * iso_dirs(n, rng)))


def totals(t):
    tt = {s: t.totals(s) for s in ("ph", "e-", "e+")}
    E = sum(v["E"] for v in tt.values()); P = np.sum([v["P"] for v in tt.values()], 0)
    return tt, E, P


def boost(beta, p0, p):
    b2 = beta @ beta; g = 1 / np.sqrt(1 - b2); u = g * beta
    bp = p @ beta; fac = (g - 1) / b2 if b2 > 1e-24 else 0.5
    return g * p0 - p @ u, p - np.outer(p0, u) + np.outer(fac * bp, beta)


def mean_gamma(T):
    return 3 * T + kv(1, 1 / T) / kv(2, 1 / T)


def equilibrium(E, N, L=0.0):
    """closed box, 2->2 channels: T, W_ph, W_-, W_+ from detailed balance of gg <-> e-e+"""
    def resid(T):
        q = 0.25 * (kv(2, 1 / T) / T ** 2) ** 2
        if abs(L) < 1e-12:
            Wg = N / (1 + 2 * np.sqrt(q)); Wm = Wp = np.sqrt(q) * Wg
        else:
            a = 1 - 4 * q; b = L + 4 * q * (N - L); c = -q * (N - L) ** 2
            Wp = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a); Wm = Wp + L; Wg = N - Wp - Wm
        return Wg * 3 * T + (Wm + Wp) * mean_gamma(T) - E, (Wg, Wm, Wp)
    T = brentq(lambda T: resid(T)[0], 0.02, 50.0)
    return (T,) + resid(T)[1]


def A_minus_1(theta):
    """mean amplification per scattering - 1, isotropic MJ electrons, Thomson limit"""
    return 16 * theta * theta + 4 * theta * kv(1, 1 / theta) / kv(2, 1 / theta)


def rate_ann(um, up):
    gm = np.sqrt(1 + np.sum(um * um, -1)); gp = np.sqrt(1 + np.sum(up * up, -1))
    qe = np.maximum(gm * gp - np.sum(um * up, -1), 1 + 1e-14)
    g2 = (qe + 1) / 2; b = np.sqrt(1 - 1 / g2)
    Lb = np.where(b > 1e-3, 2 * np.arctanh(np.minimum(b, 1 - 1e-16)) / np.maximum(b, 1e-300), 2 + 2 * b * b / 3)
    return (3 / 32) / (gm * gp) * ((2 + 2 / g2 - 1 / g2 ** 2) * Lb - 2 - 2 / g2)


def mj_ann_average(theta, rng, n=400000):
    return float(np.mean(2 * rate_ann(mj_momenta(n, theta, rng), mj_momenta(n, theta, rng))))


def sigma_gg(q):
    if q <= 2 * (1 + 1e-12):
        return 0.0
    b = np.sqrt(1 - 2 / q); Lb = 2 * np.arctanh(min(b, 1 - 1e-16))
    return (3 / 16) * (1 - b * b) * ((3 - b ** 4) * Lb - 2 * b * (2 - b * b))


def rho_gg(x1, x2):
    """angle-averaged sigma-hat F for isotropic monoenergetic photons"""
    f = lambda mu: sigma_gg(x1 * x2 * (1 - mu)) * (1 - mu)
    mu_thr = 1 - 2 / (x1 * x2)
    if mu_thr >= 1:
        return 0.0
    return 0.5 * quad(f, -1, min(mu_thr, 1), limit=200)[0]


# ======================================================================================= tests
class qed_twobody(CheckMixin, unittest.TestCase):

    def test8_rate_exactness(self):
        """first-step event counts against the exact pair sum, for several steps"""
        Nph, Npm, nseed = 8000, 2000, 60  # the seed mean is good to 0.7 %
        for dt in (1e-2, 5e-3, 2e-3):
            base = make_tile(1, dt=dt); init_pairs(base, Nph, Npm, 0.05, 2.0, np.random.default_rng(1))
            exp = sum(base.pair_sum(c) for c in ("compton-", "compton+", "pairann", "photann"))
            meas = []
            for sd in range(nseed):
                t = make_tile(1 + 100 * sd, dt=dt); init_pairs(t, Nph, Npm, 0.05, 2.0, np.random.default_rng(1))
                t.interact(); meas.append(sum(t.counts().values()))
            meas = np.array(meas)
            self.check(f"events/expected dt={dt:.0e} (+- {meas.std()/np.sqrt(nseed)/exp:.3f})", 1.0, meas.mean() / exp, 0.03)

    def test1_boosted_box(self):
        """closed box boosted along z: the equilibrium follows from the invariant mass"""
        beta, steps, Nph, Npm = 0.5, 3000, 8000, 2000
        dt = 5e-3
        rng = np.random.default_rng(1)
        t = make_tile(1, dt=dt)
        bvec = np.array([0.0, 0.0, beta]); G = 1 / np.sqrt(1 - beta * beta)
        u, w = wien_pop(Nph, 0.05, rng); _, u = boost(-bvec, np.linalg.norm(u, axis=1), u); set_pop(t, "ph", u, w)
        for s in ("e-", "e+"):
            u, w = mj_pop(Npm, 2.0, rng); _, u = boost(-bvec, np.sqrt(1 + np.sum(u * u, 1)), u); set_pop(t, s, u, w)
        tt, E0, P0 = totals(t)
        N0 = sum(v["W"] for v in tt.values()); L0 = tt["e-"]["W"] - tt["e+"]["W"]
        M = np.sqrt(E0 * E0 - P0 @ P0); T, Wg, Wm, Wp = equilibrium(M, N0, L0)
        rec = dict(Wg=[], Wpm=[], dP=[], dE=[], x=[], g=[])
        for n in range(steps + 1):
            tt, E, P = totals(t)
            rec["Wg"].append(tt["ph"]["W"]); rec["Wpm"].append(tt["e-"]["W"] + tt["e+"]["W"])
            rec["dP"].append(np.linalg.norm(P - P0) / E0); rec["dE"].append(abs(E - E0) / E0)
            rec["x"].append(tt["ph"]["E"] / tt["ph"]["W"]); rec["g"].append((tt["e-"]["E"] + tt["e+"]["E"]) / (tt["e-"]["W"] + tt["e+"]["W"]))
            if n < steps:
                t.interact()
        rec = {k: np.array(v) for k, v in rec.items()}
        late = slice(int(0.75 * steps), None)
        self.check("W_gamma", Wg, rec["Wg"][late].mean(), 0.02)
        self.check("W_pm", Wm + Wp, rec["Wpm"][late].mean(), 0.02)
        self.check("<x>_lab", 3 * T * G, rec["x"][late].mean(), 0.02)
        self.check("<gamma>_lab", G * mean_gamma(T), rec["g"][late].mean(), 0.02)
        self.check("|dP|/E0", 0.0, rec["dP"].max(), 0.03)
        self.check("|dE|/E0", 0.0, rec["dE"].max(), 0.03)

    def compton_box(self, seed, uph, wph, Ne, theta, dt=5e-3):
        t = make_tile(seed, dt=dt, channels=["compton-"])
        rng = np.random.default_rng(seed + 3)
        set_pop(t, "ph", uph, wph); u, w = mj_pop(Ne, theta, rng); set_pop(t, "e-", u, w)
        set_pop(t, "e+", np.zeros((0, 3)), np.zeros(0))
        return t

    def test2_compton_thermostat(self):
        dt = 5e-3
        # (a) Wien relaxation at Theta = 0.1 from monoenergetic x0 = 0.01
        theta, x0, steps, Nph, Ne = 0.1, 0.01, 3000, 8000, 4000
        rng = np.random.default_rng(12)
        u, w = mono_pop(Nph, x0, 0.0, rng); t = self.compton_box(1, u, w, Ne, theta)
        xm = []
        for n in range(steps + 1):
            tt = t.totals("ph"); xm.append(tt["E"] / tt["W"])
            if n < steps:
                thermostat(t, "e-", theta, rng); t.interact()
        xm = np.array(xm)
        self.check("<x>/3Theta (Wien)", 1.0, xm[int(0.75 * steps):].mean() / (3 * theta), 0.03)
        # (b) Kompaneets log-normal: Theta = 0.01, x0 = 3e-5
        theta, x0, steps = 0.01, 3e-5, 600
        rng = np.random.default_rng(14)
        u, w = mono_pop(Nph, x0, 0.0, rng); t = self.compton_box(1, u, w, Ne, theta)
        tau_step = Ne * WL / N2(dt)
        for n in range(steps):
            thermostat(t, "e-", theta, rng); t.interact()
        x, ww = energies(t, "ph"); lx = np.log(x)
        mean = np.average(lx, weights=ww); var = np.average((lx - mean) ** 2, weights=ww); y = steps * tau_step * theta
        self.check("Kompaneets <ln x> drift / 3y", 1.0, (mean - np.log(x0)) / (3 * y), 0.08)
        self.check("Kompaneets Var(ln x) / 2y", 1.0, var / (2 * y), 0.10)
        # (c) amplification per scattering: photons reset to x0 every step, energy gain per moved weight
        # (the Doppler shift of one event scatters by ~2 beta, so 40 steps make the mean good to 2 %)
        x0, Nc, steps = 1e-6, 40000, 40
        for th in (0.03, 0.1, 0.3, 1.0, 3.0):
            rng = np.random.default_rng(18)
            u, w = mono_pop(Nc, x0, 0.0, rng); t = self.compton_box(1, u, w, Ne, th)
            gain = 0.0; wsc = 0.0
            for n in range(steps):
                reset_photons(t, x0, rng); thermostat(t, "e-", th, rng)
                E0 = t.totals("ph")["E"]; w0 = t.moved()["compton-"]; t.interact()
                gain += t.totals("ph")["E"] - E0; wsc += t.moved()["compton-"] - w0
            self.check(f"A-1 ratio Theta={th}", 1.0, gain / (x0 * wsc) / A_minus_1(th), 0.04)

    def test3_ic_cooling(self):
        """gamma(t) of monoenergetic electrons in a replenished soft Wien bath, against the KN-corrected law"""
        dt = 5e-3
        g0, theta, Ne, Nph, steps = 10.0, 1e-4, 2000, 20000, 300
        xmean = 3 * theta
        rng = np.random.default_rng(24)
        ue, we = mono_pop(Ne, g0, 1.0, rng); uph, wph = wien_pop(Nph, theta, rng)
        t = make_tile(1, dt=dt, channels=["compton-"]); set_pop(t, "ph", uph, wph); set_pop(t, "e-", ue, we)
        set_pop(t, "e+", np.zeros((0, 3)), np.zeros(0))
        tau_step = wph.sum() / N2(dt); gm = []
        for n in range(steps + 1):
            e, w = energies(t, "e-"); gm.append(np.average(e, weights=w))
            if n < steps:
                t.interact()
                reset_photons(t, rng.gamma(3.0, theta, int(t.totals("ph")["N"])), rng)
        gm = np.array(gm); tau = np.arange(steps + 1) * tau_step
        sol = solve_ivp(lambda tt, g: [-(4 / 3) * (g[0] ** 2 - 1) * xmean * (1 - 6.3 * g[0] * 4 * theta)], [0, tau[-1]], [g0], t_eval=tau, rtol=1e-8)
        self.check("gamma(end)/gamma_KN", 1.0, gm[-1] / sol.y[0][-1], 0.02)

    def test4_annihilation(self):
        dt = 5e-3
        Npm, steps = 4000, 200
        # (a) cold decay: W(n) = W_0/(1 + (3/8) W_0 n/N_2)
        rng = np.random.default_rng(32)
        t = make_tile(1, dt=dt, channels=["pairann"])
        um, wm = mj_pop(Npm, 1e-3, rng); up, wp = mj_pop(Npm, 1e-3, rng)
        set_pop(t, "e-", um, wm); set_pop(t, "e+", up, wp); set_pop(t, "ph", np.zeros((0, 3)), np.zeros(0))
        W = []
        for n in range(steps + 1):
            W.append(t.totals("e-")["W"])
            if n < steps:
                t.interact()
        W = np.array(W); W_th = W[0] / (1 + (3 / 8) * (W[0] / N2(dt)) * np.arange(steps + 1))
        self.check("cold W(end)/W_th", 1.0, W[-1] / W_th[-1], 0.05)
        # the annihilation line at Theta = 0.01
        rng = np.random.default_rng(38); th = 0.01
        t = make_tile(1, dt=dt, channels=["pairann"])
        um, wm = mj_pop(Npm, th, rng); up, wp = mj_pop(Npm, th, rng)
        set_pop(t, "e-", um, wm); set_pop(t, "e+", up, wp); set_pop(t, "ph", np.zeros((0, 3)), np.zeros(0))
        for n in range(60):
            t.interact()
        x, w = energies(t, "ph"); xc = 1 + 1.5 * th
        sig = np.sqrt(np.average((x - xc) ** 2, weights=w))
        self.check("line width / sqrt(Theta/2)", 1.0, sig / np.sqrt(th / 2), 0.05)
        # (b) thermal rates: leptons thermostatted while they deplete
        for th in (0.05, 0.2, 1.0):
            rng = np.random.default_rng(42)
            t = make_tile(1, dt=dt, channels=["pairann"])
            um, wm = mj_pop(2 * Npm, th, rng); up, wp = mj_pop(2 * Npm, th, rng)
            set_pop(t, "e-", um, wm); set_pop(t, "e+", up, wp); set_pop(t, "ph", np.zeros((0, 3)), np.zeros(0))
            denom = 0.0
            for n in range(600):  # 2000-3000 events: the Poisson scatter of the rate is 2 %
                denom += t.totals("e-")["W"] * t.totals("e+")["W"] / N2(dt)
                thermostat(t, "e-", th, rng); thermostat(t, "e+", th, rng); t.interact()
            rho = t.counts()["pairann"] * WL / denom
            self.check(f"thermal rate ratio Theta={th}", 1.0, rho / mj_ann_average(th, np.random.default_rng(5)), 0.04)

    def test5_photon_photon(self):
        dt = 5e-3
        x1, x2, n1, n2, steps = 0.6, 2.5, 6000, 6000, 300
        r12, r22 = rho_gg(x1, x2), rho_gg(x2, x2)
        rng = np.random.default_rng(44)
        t = make_tile(1, dt=dt, channels=["photann"])
        u1, w1 = mono_pop(n1, x1, 0.0, rng, w=WL); u2, w2 = mono_pop(n2, x2, 0.0, rng, w=WL)
        u1 *= (1 + 1e-3 * (rng.random(n1) - 0.5))[:, None]; u2 *= (1 + 1e-3 * (rng.random(n2) - 0.5))[:, None]
        set_pop(t, "ph", np.concatenate([u1, u2]), np.concatenate([w1, w2]))
        set_pop(t, "e-", np.zeros((0, 3)), np.zeros(0)); set_pop(t, "e+", np.zeros((0, 3)), np.zeros(0))
        W10, W20 = w1.sum(), w2.sum()
        def ode(tt, y):
            W1, W2, Wp = y; c12 = W1 * W2 * r12 / N2(dt); c22 = W2 * W2 * r22 / N2(dt)
            return [-c12, -c12 - c22, c12 + 0.5 * c22]
        sol = solve_ivp(ode, [0, steps], [W10, W20, 0.0], t_eval=[steps], rtol=1e-8)
        for n in range(steps):
            t.interact()
        # every pair is born at one position: electrons and positrons coincide, so no charge
        # density is added and div E is what the fields had
        pm = np.sort(np.stack(t.get_positions(t.species("e-")), 1), axis=0)
        pp = np.sort(np.stack(t.get_positions(t.species("e+")), 1), axis=0)
        self.assertGreater(len(pm), 100)
        np.testing.assert_array_equal(pm, pp)
        x, w = energies(t, "ph")
        # the weight the reactions moved, which is what the ~3300 events make good to 2 %
        self.check("W1 loss/ODE", 1.0, (W10 - w[x < 1.0].sum()) / (W10 - sol.y[0][-1]), 0.05)
        self.check("W2 loss/ODE", 1.0, (W20 - w[x >= 1.0].sum()) / (W20 - sol.y[1][-1]), 0.05)
        self.check("W_pm(end)/ODE", 1.0, t.totals("e-")["W"] / sol.y[2][-1], 0.05)


if __name__ == "__main__":
    unittest.main()
