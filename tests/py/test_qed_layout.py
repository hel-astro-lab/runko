# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import numpy as np
import runko

NC = 8
CFL = 0.5
NBIN = 8
NRC = 2
SLACK = 0.5
SMIN = 4
WL = 1e-3


def make_tile(positron_first=False):
    config = runko.Configuration(None)
    config.n_tiles = [1, 1, 1]
    config.n_cells_per_tile = [NC, NC, NC]
    config.cfl = CFL
    config.field_propagator = "fdtd2"
    if positron_first:
        config.q0, config.m0, config.w0 = 1, 1, WL
        config.q1, config.m1, config.w1 = -1, 1, WL
    else:
        config.q0, config.m0, config.w0 = -1, 1, WL
        config.q1, config.m1, config.w1 = 1, 1, WL
    config.q2, config.m2 = 0, 0
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    config.qed_H_lamC = 2.33e17
    config.qed_nbin = NBIN
    config.qed_n_reaction_cells = NRC
    config.qed_slack = SLACK
    config.qed_smin = SMIN
    config.qed_seed = 3
    tile = runko.qed.threeD.Tile((0, 0, 0), config)
    zero = lambda x, y, z: (0, 0, 0)
    tile.set_EBJ(zero, zero, zero)
    return config, tile


def populate(tile, rng, nph=3000, nlep=1000):
    x = np.exp(rng.uniform(np.log(1e-4), np.log(10.0), nph))
    d = rng.normal(size=(nph, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    tile.set_particles("ph", [tuple(v) for v in x[:, None] * d], list(rng.uniform(0.3, 3.0, nph) * WL))
    for s in ("e-", "e+"):
        p = np.exp(rng.uniform(np.log(0.05), np.log(30.0), nlep))
        d = rng.normal(size=(nlep, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
        tile.set_particles(s, [tuple(v) for v in p[:, None] * d], [WL] * nlep)


def energy(tile, s):
    m = 0.0 if s == "ph" else 1.0
    vx, vy, vz = tile.get_velocities(tile.species(s))
    return np.sqrt(m * m + np.asarray(vx) ** 2 + np.asarray(vy) ** 2 + np.asarray(vz) ** 2)


def check_layout(test, tile, s, mins, maxs):
    """every live LP sits in the key of its cell and energy bin, inside the bounds of the key"""
    sp = tile.species(s)
    keyed = s == "ph"
    m = 0.0 if s == "ph" else 1.0
    keys = np.asarray(tile.get_slot_keys(s))
    tb = tile.get_tables(s)
    e = energy(tile, s)
    w = np.asarray(tile.get_weights(sp))
    x, y, z = (np.asarray(a) for a in tile.get_positions(sp))
    n = len(e)
    test.assertEqual(len(keys), n)
    test.assertTrue(np.all(keys >= 0))
    nkeys = (2 if keyed else 1) * NRC ** 3 * NBIN
    test.assertEqual(len(tb["b"]), nkeys)
    # counts
    cnt = np.bincount(keys, minlength=nkeys)
    np.testing.assert_array_equal(cnt, tb["fill"])
    test.assertTrue(np.all(tb["cap"] >= tb["fill"] + np.maximum(np.ceil(SLACK * tb["fill"]), SMIN)))
    np.testing.assert_array_equal(tb["b"], np.concatenate([[0], np.cumsum(tb["cap"])[:-1]]))
    # bounds
    test.assertTrue(np.all(e <= tb["emax"][keys] * (1 + 1e-6)))
    test.assertTrue(np.all(e >= tb["emin"][keys] * (1 - 1e-6)))
    test.assertTrue(np.all(w <= tb["wmax"][keys] * (1 + 1e-6)))
    live = tb["fill"] > 0
    test.assertTrue(np.all(tb["emax"][live] > 0))
    test.assertTrue(np.all(np.isinf(tb["emin"][~live])))
    # cell and bin of the key
    kk = keys // 2 if keyed else keys
    cell = kk // NBIN
    ebin = kk % NBIN
    if NBIN > 1:
        v = np.maximum(e - m, 1e-30)
        b_exp = np.clip(np.floor((np.log(v) - tb["llo"]) / tb["dl"]).astype(int), 0, NBIN - 1)
        np.testing.assert_array_equal(ebin, b_exp)
    ci = np.clip(np.floor((x - mins[0]) / (maxs[0] - mins[0]) * NRC).astype(int), 0, NRC - 1)
    cj = np.clip(np.floor((y - mins[1]) / (maxs[1] - mins[1]) * NRC).astype(int), 0, NRC - 1)
    ck = np.clip(np.floor((z - mins[2]) / (maxs[2] - mins[2]) * NRC).astype(int), 0, NRC - 1)
    np.testing.assert_array_equal(cell, (ci * NRC + cj) * NRC + ck)
    # the slots are grouped: keys nondecreasing in slot order
    test.assertTrue(np.all(np.diff(keys) >= 0))
    # photon bits both occur
    if keyed:
        test.assertTrue(0 < np.count_nonzero(keys % 2) < n)


class qed_layout(unittest.TestCase):

    def test_box_and_roles(self):
        for pf in (False, True):
            config, tile = make_tile(pf)
            b = tile.box()
            self.assertAlmostEqual(b["dt"], CFL / NC)
            self.assertAlmostEqual(b["N2"], NC / CFL)
            self.assertAlmostEqual(b["N1"], 1.0 / (7.2973525693e-3 * b["dt"] * config.qed_H_lamC), places=6)
            self.assertEqual(b["w_lep"], WL)
            self.assertEqual(b["b_perp"], 0.0)
            self.assertEqual(tile.species("ph"), 2)
            self.assertEqual(tile.species("e-"), 1 if pf else 0)
            self.assertEqual(tile.species("e+"), 0 if pf else 1)
            self.assertFalse(tile.laid())

    def test_layout_tables(self):
        config, tile = make_tile()
        rng = np.random.default_rng(5)
        populate(tile, rng)
        t0 = {s: tile.totals(s) for s in ("ph", "e-", "e+")}
        tile.sort_particles()
        self.assertTrue(tile.laid())
        mins, maxs = (0.0, 0.0, 0.0), (NC, NC, NC)
        for s in ("ph", "e-", "e+"):
            check_layout(self, tile, s, mins, maxs)
            t1 = tile.totals(s)
            self.assertEqual(t1["N"], t0[s]["N"])
            self.assertAlmostEqual(t1["E"], t0[s]["E"], delta=1e-6 * t0[s]["E"])
            self.assertAlmostEqual(t1["W"], t0[s]["W"], delta=1e-6 * t0[s]["W"])
        self.assertAlmostEqual(tile.tau_T(), 2000 * WL, places=9)

    def test_relayout_is_consistent(self):
        config, tile = make_tile()
        rng = np.random.default_rng(6)
        populate(tile, rng)
        tile.sort_particles()
        f0 = {s: tile.get_tables(s)["fill"].copy() for s in ("ph", "e-", "e+")}
        n0 = {s: tile.totals(s)["N"] for s in ("ph", "e-", "e+")}
        tile.sort_particles()
        mins, maxs = (0.0, 0.0, 0.0), (NC, NC, NC)
        for s in ("ph", "e-", "e+"):
            check_layout(self, tile, s, mins, maxs)
            self.assertEqual(tile.totals(s)["N"], n0[s])
            if s != "ph":  # the photon bit is redrawn every epoch
                np.testing.assert_array_equal(tile.get_tables(s)["fill"], f0[s])

    def test_tail_after_injection(self):
        config, tile = make_tile()
        rng = np.random.default_rng(7)
        populate(tile, rng)
        tile.sort_particles()
        n0 = tile.totals("ph")["N"]
        # particles appended outside the layout are inert until the next sort
        tile.inject(tile.species("ph"), [runko.pic.threeD.ParticleState(pos=(1.0, 2.0, 3.0), vel=(0.5, 0, 0), wgt=1e-3)
                                          for _ in range(50)])
        self.assertEqual(tile.totals("ph")["N"], n0 + 50)
        tile.sort_particles()
        check_layout(self, tile, "ph", (0.0, 0.0, 0.0), (NC, NC, NC))
        self.assertEqual(tile.totals("ph")["N"], n0 + 50)

    def test_empty_species(self):
        config, tile = make_tile()
        rng = np.random.default_rng(8)
        populate(tile, rng, nph=0, nlep=10)
        tile.sort_particles()
        tb = tile.get_tables("ph")
        self.assertTrue(np.all(tb["fill"] == 0))
        self.assertTrue(np.all(tb["cap"] == SMIN))
        self.assertEqual(len(tile.get_slot_keys("ph")), 0)
        check_layout(self, tile, "e-", (0.0, 0.0, 0.0), (NC, NC, NC))


if __name__ == "__main__":
    unittest.main()
