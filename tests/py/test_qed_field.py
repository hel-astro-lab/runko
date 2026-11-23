# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""The QED channels on the local lattice field: every particle's quantum parameter comes from
the interpolated E, B at its position, scaled so that the reference field qed_B0 (code units)
has the quantum parameter of the magnetic compactness qed_ell_B."""

import unittest
import numpy as np
import runko
from tests.qed_test_helpers import H, WL, NC, make_tile, iso_dirs, set_pop, empty, state, CheckMixin

ph = runko.qed.physics


def make(seed, N1, ell_B, B0, channels, Bz_of_x, acc_syn=0.0):
    """a tile whose lattice field is Bz(x) (code units), the channels reading it"""
    dt = 1.0 / (ph.ALPHA * N1 * H)
    return make_tile(seed, dt, Bz=Bz_of_x, ell_B=ell_B, channels=channels, acc_syn=acc_syn, qed_B0=B0), dt


def positions(t, s):
    return np.stack([np.asarray(a) for a in t.get_positions(t.species(s))], 1)


class qed_field(CheckMixin, unittest.TestCase):

    def test_emission_follows_the_local_field(self):
        """gamma = 30 leptons, isotropic, in a lattice whose Bz is B0 for x < NC/2 and B0/4 beyond:
        each lepton emits at the rate of its own b_eff = b0 |beta x B|/B0, so the emission
        counted per half agrees with the per-lepton expectation of that half"""
        g0, N1, B0, ell_B = 30.0, 0.02, 2.0, 10.0
        b0 = float(ph.b_from_ellB(ell_B, H))
        rng = np.random.default_rng(3)
        n = 20000
        u = np.sqrt(g0 * g0 - 1) * iso_dirs(n, rng)
        nseed, ev = 6, np.zeros((2, 6))
        for sd in range(nseed):
            t, dt = make(1 + 10 * sd, N1, ell_B, B0, ["syn-"], lambda x: np.where(x < NC / 2, B0, B0 / 4))
            set_pop(t, "e-", u, WL * np.ones(n)); empty(t, ["e+", "ph"])
            pos = positions(t, "e-")
            uu, gam, _ = state(t, "e-")
            # the expectation per lepton: b_eff = b0 |beta x z| B(x)/B0 = b0 sin(pitch) B(x)/B0
            beta = uu / gam[:, None]
            sina = np.sqrt(beta[:, 0] ** 2 + beta[:, 1] ** 2) / np.linalg.norm(beta, axis=1)
            Bx = np.where(pos[:, 0] < NC / 2, B0, B0 / 4)
            b_eff = b0 * sina * Bx / B0
            lam = ph.rate_syn(uu[:, 0], uu[:, 1], uu[:, 2], b_eff, 0.0, False) / N1
            p = 1 - np.exp(-lam)
            left = pos[:, 0] < NC / 2
            exp = np.array([p[left].sum(), p[~left].sum()])
            t.interact()
            # the emitted photons sit at their parents' positions: count them per half
            px = positions(t, "ph")[:, 0]
            ev[:, sd] = [np.count_nonzero(px < NC / 2), np.count_nonzero(px >= NC / 2)]
        for half, name in ((0, "strong half"), (1, "weak half")):
            m = ev[half].mean(); e = ev[half].std() / np.sqrt(nseed)
            self.check(f"emissions in the {name} / expected (+- {e/exp[half]:.3f})", 1.0, m / exp[half], max(0.03, 4 * e / exp[half]))
        self.assertGreater(ev[0].mean(), 2 * ev[1].mean())

    def test_pair_production_follows_the_local_field(self):
        """x = 20 photons along z in a lattice Bz(x): no conversions where B is along k,
        the expected rate where the field is transverse"""
        N1, B0, ell_B, x0 = 0.5, 1.0, 2.4323, 20.0
        H = 1e5
        dt = 1.0 / (ph.ALPHA * N1 * H)
        b0 = float(ph.b_from_ellB(ell_B, H))
        rng = np.random.default_rng(5)
        n = 20000
        c = runko.Configuration(None)
        c.n_tiles = [1, 1, 1]; c.n_cells_per_tile = [NC, NC, NC]; c.cfl = dt * NC; c.field_propagator = "fdtd2"
        c.q0, c.m0, c.w0 = -1, 1, WL; c.q1, c.m1, c.w1 = 1, 1, WL; c.q2, c.m2 = 0, 0
        c.particle_pusher = ["boris", "boris", "ballistic"]; c.field_interpolator = "linear_1st"; c.current_depositer = "zigzag_1st"
        c.qed_H_lamC = H; c.qed_ell_B = ell_B; c.qed_B0 = B0
        c.qed_nbin, c.qed_nsort, c.qed_slack, c.qed_smin = 8, 1, 0.1, 4; c.qed_a_res = 0.0; c.qed_seed = 9
        c.qed_channels = ["bw"]
        t = runko.qed.threeD.Tile((0, 0, 0), c)
        zero = lambda x, y, z: (0, 0, 0)
        # the field is along z on the left half, along x (transverse to k = x0 z) on the right
        t.set_EBJ(zero, lambda x, y, z: (np.where(x < NC / 2, 0.0, B0), 0, np.where(x < NC / 2, B0, 0.0)), zero)
        # photons two cells away from the half-boundary, where the interpolated field is pure
        xs = np.where(rng.random(n) < 0.5, rng.uniform(2.0, NC / 2 - 2.0, n), rng.uniform(NC / 2 + 2.0, NC - 2.0, n))
        t.inject(t.species("ph"), [runko.pic.threeD.ParticleState(pos=(x, rng.uniform(1, NC - 1), rng.uniform(1, NC - 1)),
                                                                    vel=(0.0, 0.0, x0), wgt=WL) for x in xs])
        empty(t, ["e-", "e+"])
        right = xs >= NC / 2
        lam = float(ph.rate_bw(0.0, 0.0, x0, b0)) / N1
        exp = (1 - np.exp(-lam)) * np.count_nonzero(right)
        t.interact()
        pe = positions(t, "e-")[:, 0]
        self.assertEqual(np.count_nonzero(pe < NC / 2), 0, "no conversion with the field along k")
        nconv = np.count_nonzero(pe >= NC / 2)
        self.check("conversions in the transverse half / expected", 1.0, nconv / exp, max(0.05, 4 / np.sqrt(exp)))

    def test_zero_field_is_silent(self):
        t, dt = make(2, 0.02, 10.0, 2.0, ["syn-", "ssa-", "bw"], lambda x: 0.0 * x)
        rng = np.random.default_rng(6)
        n = 5000
        set_pop(t, "e-", np.sqrt(30.0 ** 2 - 1) * iso_dirs(n, rng), WL * np.ones(n))
        set_pop(t, "ph", 20.0 * iso_dirs(n, rng), WL * np.ones(n)); empty(t, ["e+"])
        for _ in range(3):
            t.interact()
        self.assertEqual(sum(t.counts().values()), 0.0)


if __name__ == "__main__":
    unittest.main()
