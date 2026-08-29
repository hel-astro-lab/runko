# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""What the qed tests share: the closed-box tile from its physical inputs, the population
setters and samplers, the state readers, and the relative-error check."""

import numpy as np
import runko

ph = runko.qed.physics

WL = 1e-3          # lepton LP weight, Thomson depths
H = 4.111e9        # tile edge in Compton wavelengths: ell_B = 10 -> b_perp = 1e-3
NC = 64            # cells per tile edge; dt = cfl/NC crossings


def make_conf(seed, dt, H_lamC=H, w_lep=WL, ell_B=0.0, x_min=0.0, channels=None, nbin=8, nsort=1, slack=0.1,
              a_res=0.0, x_ref=1e-2, kappa=0.1, acc_syn=3.0, smin=4, nsub_max=8, nc=NC, **extra):
    """one cubic tile of nc^3 cells, electrons (w0), positrons (w1) and photons, stepping dt crossings;
    extra keys go to the config as they are"""
    c = runko.Configuration(None)
    c.n_tiles = [1, 1, 1]
    c.n_cells_per_tile = [nc, nc, nc]
    c.cfl = dt * nc
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
    for k, v in extra.items():
        setattr(c, k, v)
    return c


def make_tile(seed, dt, Bz=None, **kw):
    """the tile of make_conf with a zero lattice, or with Bz(x) (code units) along z"""
    t = runko.qed.threeD.Tile((0, 0, 0), make_conf(seed, dt, **kw))
    zero = lambda x, y, z: (0, 0, 0)
    t.set_EBJ(zero, zero if Bz is None else (lambda x, y, z: (0, 0, Bz(x))), zero)
    return t


def box_N1(seed, N1, ell_B=10.0, H_lamC=H, x_min=0.0, **kw):
    """the box whose one-body clock is N1: the step follows from the size"""
    return make_tile(seed, 1.0 / (ph.ALPHA * N1 * H_lamC), H_lamC=H_lamC, ell_B=ell_B, x_min=x_min, **kw)


def rows(u):
    return [tuple(map(float, v)) for v in np.asarray(u).reshape(-1, 3)]


def set_pop(t, s, u, w):
    t.set_particles(s, rows(u), [float(x) for x in np.asarray(w)])


def empty(t, species):
    for s in species:
        t.set_particles(s, [], [])


def iso_dirs(n, rng):
    cz = -1 + 2 * rng.random(n); phi = 2 * np.pi * rng.random(n); sz = np.sqrt(1 - cz * cz)
    return np.stack([sz * np.cos(phi), sz * np.sin(phi), cz], 1)


def mj_momenta(n, theta, rng):
    """n isotropic Maxwell-Juettner momenta at temperature theta, as an (n, 3) array"""
    return np.stack(runko.sample_boosted_juttner_synge(n, theta, beta=0, gen=rng), 1)


def spread(n, rng):
    """photon weights spread over half a decade about the lepton weight"""
    return WL * np.exp(rng.uniform(np.log(0.3), np.log(3.0), n))


def mono_pop(n, e, mass, rng, w=None):
    """n isotropic LPs of energy e; leptons at the lepton weight, photons at spread weights, or all at w"""
    p = np.sqrt(max(e * e - mass * mass, 0.0)) if mass > 0 else e
    ww = (WL * np.ones(n) if mass > 0 else spread(n, rng)) if w is None else w * np.ones(n)
    return p * iso_dirs(n, rng), ww


def state(t, s):
    """momenta (n, 3), energies and weights of the live LPs of a species"""
    m = 0.0 if s == "ph" else 1.0
    vx, vy, vz = (np.asarray(a) for a in t.get_velocities(t.species(s)))
    u = np.stack([vx, vy, vz], 1)
    return u, np.sqrt(m * m + np.sum(u * u, 1)), np.asarray(t.get_weights(t.species(s)))


def totals(t):
    """the totals of every species, and the summed energy"""
    tt = {s: t.totals(s) for s in ("ph", "e-", "e+")}
    return tt, sum(v["E"] for v in tt.values())


class CheckMixin:
    """relative-error assertion of a measured value against its prediction"""

    def check(self, name, pred, meas, tol):
        err = abs(meas - pred) / (abs(pred) if pred != 0 else 1.0)
        self.assertLessEqual(err, tol, f"{name}: predicted {pred:.5g}, measured {meas:.5g}, tol {tol}")
