# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import re
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


def make_tile(positron_first=False, proton=False, bad_mass=None, idx=(0, 0, 0), **extra):
    config = runko.Configuration(None)
    config.n_tiles = [2, 1, 1]
    config.n_cells_per_tile = [NC, NC, NC]
    config.cfl = CFL
    config.field_propagator = "fdtd2"

    if proton:
        # an ion species first: it must be skipped, not mistaken for a lepton
        config.q0, config.m0 = 1, 1836
        config.q1, config.m1, config.w1 = 1, 1, WL
        config.q2, config.m2, config.w2 = -1, 1, WL
        config.q3, config.m3, config.n_payload3 = 0, 0, 1
    elif positron_first:
        config.q0, config.m0, config.w0 = 1, 1, WL
        config.q1, config.m1, config.w1 = -1, 1, WL
        config.q2, config.m2, config.n_payload2 = 0, 0, 1
    else:
        config.q0, config.m0, config.w0 = -1, 1, WL
        config.q1, config.m1, config.w1 = 1, 1, WL
        config.q2, config.m2, config.n_payload2 = 0, 0, 1
    if bad_mass is not None:
        config.q0, config.m0 = -1, bad_mass

    config.particle_pusher = ["boris"] * (3 if proton else 2) + ["ballistic"]
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"

    config.qed_H_lamC = 2.33e17
    config.qed_nbin = NBIN
    config.qed_n_reaction_cells = NRC
    config.qed_slack = SLACK
    config.qed_smin = SMIN
    config.qed_seed = 3
    for k, v in extra.items():  # further config keys as they are
        setattr(config, k, v)

    tile = runko.qed.threeD.Tile(idx, config)
    zero = lambda x, y, z: (0, 0, 0)
    tile.set_EBJ(zero, zero, zero)

    return config, tile


def populate(tile, rng, nph=3000, nlep=1000):
    x = np.exp(rng.uniform(np.log(1e-4), np.log(10.0), nph))
    d = rng.normal(size=(nph, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    tile.set_particles("ph", [tuple(v) for v in x[:, None] * d],
                       list(rng.uniform(0.3, 3.0, nph) * WL))

    for s in ("e-", "e+"):
        p = np.exp(rng.uniform(np.log(0.05), np.log(30.0), nlep))
        d = rng.normal(size=(nlep, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
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
    test.assertTrue(np.all(tb["cap"]
                           >= tb["fill"] + np.maximum(np.ceil(SLACK * tb["fill"]), SMIN)))
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
            self.assertAlmostEqual(b["N1"], 1.0 / (7.2973525693e-3 * b["dt"] * config.qed_H_lamC),
                                   places=6)
            self.assertEqual(b["w_lep"], WL)
            self.assertEqual(b["b_perp"], 0.0)
            self.assertEqual(tile.species("ph"), 2)
            self.assertEqual(tile.species("e-"), 1 if pf else 0)
            self.assertEqual(tile.species("e+"), 0 if pf else 1)
            self.assertFalse(tile.laid())

        # a proton species (m > 1) is skipped, whatever its charge and position in the list
        config, tile = make_tile(proton=True)
        self.assertEqual(tile.species("ph"), 3)
        self.assertEqual(tile.species("e+"), 1)
        self.assertEqual(tile.species("e-"), 2)

        # a charged species of fractional mass is neither a lepton nor an ion
        with self.assertRaises(RuntimeError):
            make_tile(bad_mass=0.5)

        # the photon weight is payload component 0; the leptons carry the species weight
        for payloads in ({"n_payload2": 0}, {"n_payload0": 1}):
            with self.assertRaises(RuntimeError):
                make_tile(**payloads)

    def test_layout_tables(self):
        config, tile = make_tile()
        rng = np.random.default_rng(5)
        populate(tile, rng)
        t0 = {s: tile.totals(s) for s in ("ph", "e-", "e+")}

        tile.layout()
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

        tile.layout()
        f0 = {s: tile.get_tables(s)["fill"].copy() for s in ("ph", "e-", "e+")}
        n0 = {s: tile.totals(s)["N"] for s in ("ph", "e-", "e+")}
        tile.layout()

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
        tile.layout()
        n0 = int(tile.totals("ph")["N"])

        tile.inject(tile.species("ph"),
                    [runko.pic.threeD.ParticleState(pos=(1.0, 2.0, 3.0), vel=(0.5, 0, 0))
                     for _ in range(50)],
                    payload=[1e-3] * 50)
        self.assertEqual(tile.totals("ph")["N"], n0 + 50)

        # pic's append writes behind the last alive slot, into the slack of the keys there and
        # past them: the tables follow the container, and the injected LPs are no members of a key
        self.assertTrue(tile.laid())
        keys = np.asarray(tile.get_slot_keys("ph"))
        tb = tile.get_tables("ph")
        self.assertEqual(len(keys), n0 + 50)
        np.testing.assert_array_equal(np.bincount(keys[:n0], minlength=len(tb["fill"])),
                                      tb["fill"])

        tile.layout()
        self.assertTrue(tile.laid())
        check_layout(self, tile, "ph", (0.0, 0.0, 0.0), (NC, NC, NC))
        self.assertEqual(tile.totals("ph")["N"], n0 + 50)

    def test_arrivals_keep_the_layout(self):
        """
        photons crossing from the neighbour tile are appended behind the last alive slot without
        invalidating the layout; the next sort keys them
        """
        import tempfile

        config = runko.Configuration(None)
        config.tile_partitioning = "catepillar_track"
        config.catepillar_track_length = 1
        config.n_laps = 1
        config.n_tiles = [2, 1, 1]
        config.n_cells_per_tile = [NC, NC, NC]
        config.cfl = 0.45
        config.field_propagator = "fdtd2"
        config.q0, config.m0, config.w0 = -1, 1, WL
        config.q1, config.m1, config.w1 = 1, 1, WL
        config.q2, config.m2, config.n_payload2 = 0, 0, 1
        config.particle_pusher = ["boris", "boris", "ballistic"]
        config.field_interpolator = "linear_1st"
        config.current_depositer = "zigzag_1st"
        config.qed_H_lamC = 2.33e17
        config.qed_nbin, config.qed_n_reaction_cells, config.qed_nsort = NBIN, NRC, 10
        config.qed_slack, config.qed_smin, config.qed_seed = SLACK, SMIN, 4
        config.io_outdir = tempfile.mkdtemp(prefix="runko-qed-arrivals-")
        tile_grid = runko.TileGrid(config)

        rng = np.random.default_rng(9)
        n_cross, n_stay = 300, 500
        tiles = {}
        for idx in tile_grid.local_tile_indices():
            tile = runko.qed.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            if idx == (0, 0, 0):
                # photons at the +x edge flying +x cross in one push; the rest stay
                w_cross = rng.uniform(0.1, 5.0, n_cross) * WL
                cross = [runko.pic.threeD.ParticleState(pos=(NC - 0.2,
                                                             1 + (NC - 2) * rng.random(),
                                                             1 + (NC - 2) * rng.random()),
                                                        vel=(1.0, 0, 0))
                         for _ in w_cross]
                w_stay = rng.uniform(0.1, 5.0, n_stay) * WL
                stay = [runko.pic.threeD.ParticleState(pos=(1 + (NC - 3) * rng.random(),
                                                            1 + (NC - 2) * rng.random(),
                                                            1 + (NC - 2) * rng.random()),
                                                       vel=tuple(rng.normal(size=3)))
                        for _ in w_stay]
                tile.inject(2, cross + stay, payload=np.concatenate([w_cross, w_stay]))
            else:
                w_stay = rng.uniform(0.1, 5.0, n_stay) * WL
                tile.inject(2,
                            [runko.pic.threeD.ParticleState(pos=(NC + 1 + (NC - 3) * rng.random(),
                                                                 1 + (NC - 2) * rng.random(),
                                                                 1 + (NC - 2) * rng.random()),
                                                            vel=tuple(rng.normal(size=3)))
                             for _ in w_stay],
                            payload=w_stay)

            tile_grid.add_tile(tile, idx)
            tiles[idx] = tile

        simulation = tile_grid.configure_simulation(config)
        # the tiles as the simulation holds them
        tiles = {tuple(tile.index): tile for tile in simulation.local_tiles()}
        for tile in simulation.local_tiles():
            tile.layout()

        def lap(x):
            x.prtcl_push()
            x.prtcl_pack_outgoing()
            x.comm_external(runko.tools.comm_mode.pic_particle)
            x.comm_local(runko.tools.comm_mode.pic_particle)

        simulation.for_each_lap(lap)

        t1 = tiles[(1, 0, 0)]
        self.assertTrue(t1.laid(), "arrivals must not invalidate the layout")
        keys = np.asarray(t1.get_slot_keys("ph"))
        self.assertEqual(len(keys), n_stay + n_cross)

        # arrivals sit behind the stayers, in slack or tail: the fill counts the stayers alone
        self.assertTrue(np.all(keys[:n_stay] >= 0))
        np.testing.assert_array_equal(np.bincount(keys[:n_stay],
                                                  minlength=len(t1.get_tables("ph")["fill"])),
                                      t1.get_tables("ph")["fill"], "arrivals are no members")
        self.assertEqual(t1.totals("ph")["N"], n_stay + n_cross)
        self.assertEqual(tiles[(0, 0, 0)].totals("ph")["N"], n_stay)

        t1.layout()
        check_layout(self, t1, "ph", (NC, 0.0, 0.0), (2 * NC, NC, NC))

    def test_arrival_in_slack_survives_products(self):
        """
        an LP appended into a key's slack keeps its slot: a product finding it there goes to the
        tail. Pair annihilation alone never touches a photon, so every photon must come out of the
        steps as it went in, the injected ones too
        """
        config, tile = make_tile(qed_channels=["pairann"], qed_nbin=1, qed_n_reaction_cells=1,
                                 qed_slack=1.0, qed_a_res=0.0, qed_nsort=1000)

        rng = np.random.default_rng(21)
        nph, nlep, k = 200, 1000, 5
        x = np.exp(rng.uniform(np.log(1e-2), 0.0, nph))
        d = rng.normal(size=(nph, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        tile.set_particles("ph", [tuple(v) for v in x[:, None] * d],
                           list(rng.uniform(0.3, 3.0, nph) * WL))
        for s in ("e-", "e+"):
            tile.set_particles(s, [tuple(v) for v in 0.3 * rng.normal(size=(nlep, 3))],
                               [WL] * nlep)

        # lays out; its products fit the slack, so the last alive photon is in the last key
        tile.interact()
        sp = tile.species("ph")

        def photons():
            u = np.stack([np.asarray(a) for a in tile.get_velocities(sp)], 1)
            return dict(zip(np.asarray(tile.get_ids(sp)),
                            zip(map(tuple, u), np.asarray(tile.get_weights(sp)))))

        before = photons()
        cap0 = tile.get_tables("ph")["cap"]
        tile.inject(sp, [runko.pic.threeD.ParticleState(pos=(1.0, 2.0, 3.0), vel=(0.5, 0, 0))
                         for _ in range(k)], payload=[7e-3] * k)
        self.assertTrue(tile.laid(), "an injection keeps the layout")

        state = photons()
        new = sorted(set(state) - set(before))
        self.assertEqual(len(new), k)

        # the injected LPs took the first slack slots of the last key, which keeps the rest
        keys = np.asarray(tile.get_slot_keys("ph"))
        np.testing.assert_array_equal(keys[-k:], 1)
        tb = tile.get_tables("ph")
        self.assertEqual(tb["fill"][1] + k, np.count_nonzero(keys == 1))
        np.testing.assert_array_equal(tb["cap"], cap0)

        held = 0.0
        for _ in range(3):
            tile.interact()
            held += float(re.search(r"\(([0-9.]+) found an arrival", tile.report()).group(1))
            ids = np.asarray(tile.get_ids(sp))
            self.assertEqual(len(ids), len(set(ids)), "an id twice")

            after = photons()
            # every photon, the injected ones too, as it was
            self.assertFalse(set(state) - set(after), "photons lost")
            self.assertEqual({i: after[i] for i in state}, state)
            state = after
        self.assertGreater(held, 0.0, "no product found an arrival in slack")

    def test_arrival_in_slack_widens_no_bound(self):
        """
        an LP appended into a key's slack is no member of the key until the fill reaches it:
        however heavy and soft, it leaves the bounds the envelopes read alone
        """
        config, tile = make_tile(qed_channels=["pairann"], qed_nbin=1, qed_n_reaction_cells=1,
                                 qed_slack=1.0, qed_a_res=0.0, qed_nsort=1000)

        rng = np.random.default_rng(22)
        nph, k = 200, 5
        d = rng.normal(size=(nph, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        tile.set_particles("ph", [tuple(v) for v in 0.1 * d], [WL] * nph)
        tile.interact()  # lays out
        sp = tile.species("ph")
        tb0 = tile.get_tables("ph")

        tile.inject(sp, [runko.pic.threeD.ParticleState(pos=(1.0, 2.0, 3.0), vel=(1e-6, 0, 0))
                         for _ in range(k)], payload=[1e3 * WL] * k)
        tile.interact()  # no layout: the step re-widens the bounds by the members

        self.assertTrue(tile.laid())
        np.testing.assert_array_equal(np.asarray(tile.get_slot_keys("ph"))[-k:], 1)
        tb1 = tile.get_tables("ph")
        for q in ("wmax", "emin", "emax", "fill"):
            np.testing.assert_array_equal(tb1[q], tb0[q])

    def test_pic_sort_on_a_qed_tile(self):
        """
        a pic sort reorders a laid-out container and drops its layout; the next step lays it out
        again (Compton only, no splitting: weight and count of every species conserved)
        """
        config, tile = make_tile(qed_channels=["compton-", "compton+"], qed_a_res=0.0,
                                 qed_n_reaction_cells=1)
        populate(tile, np.random.default_rng(13))
        tile.interact()
        self.assertTrue(tile.laid())
        t0 = {s: tile.totals(s) for s in ("ph", "e-", "e+")}

        tile.sort_particles()
        self.assertFalse(tile.laid())

        tile.interact()
        self.assertTrue(tile.laid())
        self.assertGreater(sum(tile.counts().values()), 0)
        for s in ("ph", "e-", "e+"):
            t1 = tile.totals(s)
            self.assertEqual(t1["N"], t0[s]["N"])
            self.assertAlmostEqual(t1["W"], t0[s]["W"], delta=1e-9 * t0[s]["W"])

    def test_tiles_draw_their_own_streams(self):
        """
        two tiles of one seed with the same population must not make the same random decisions:
        the layout's photon bit is drawn per pre-layout slot, which the injection order fixes
        alike in both tiles, so the bits agree only where the streams do
        """
        bits = []
        for idx in ((0, 0, 0), (1, 0, 0)):
            config, tile = make_tile(idx=idx)
            populate(tile, np.random.default_rng(12))
            tile.layout()
            ordinal = np.asarray(tile.get_ids(tile.species("ph"))) & ((1 << 40) - 1)
            bits.append(dict(zip(ordinal, np.asarray(tile.get_slot_keys("ph")) % 2)))

        common = sorted(set(bits[0]) & set(bits[1]))
        self.assertEqual(len(common), 3000)
        same = np.mean([bits[0][o] == bits[1][o] for o in common])
        self.assertLess(same, 0.6, f"independent bits agree on half, these on {same:.3f}")

    def test_species_set_together_do_not_share_positions(self):
        """set_particles draws the positions of every species from its own stream"""
        config, tile = make_tile()
        populate(tile, np.random.default_rng(10), nph=1000, nlep=1000)

        pos = {s: np.stack([np.asarray(a) for a in tile.get_positions(tile.species(s))], 1)
               for s in ("ph", "e-", "e+")}
        for a, b in (("ph", "e-"), ("e-", "e+")):
            self.assertEqual(np.count_nonzero(np.all(pos[a] == pos[b], axis=1)), 0,
                             f"{a} and {b} co-located")

    def test_qed_actions_need_a_qed_tile(self):
        """
        the qed lap actions on a plain pic tile raise instead of doing nothing; the sort still
        works
        """
        config, _ = make_tile()
        tile = runko.pic.threeD.Tile((0, 0, 0), config)

        a = runko.actions
        for sym, args in ((a.interact, ()),
                          (a.escape_photons, ()),
                          (a.heat_diffusive, (1.0,)),
                          (a.inject_pairs_plaw, (1.0, 2.0, 2.0, 10.0)),
                          (a.inject_photons, (1.0, 0.1, WL))):
            with self.assertRaises(RuntimeError):
                tile.sim_context.eval((sym, a.current_context, *args))

        tile.sort_particles()

    def test_empty_species(self):
        config, tile = make_tile()
        rng = np.random.default_rng(8)
        populate(tile, rng, nph=0, nlep=10)
        tile.layout()

        tb = tile.get_tables("ph")
        self.assertTrue(np.all(tb["fill"] == 0))
        self.assertTrue(np.all(tb["cap"] == SMIN))
        self.assertEqual(len(tile.get_slot_keys("ph")), 0)
        check_layout(self, tile, "e-", (0.0, 0.0, 0.0), (NC, NC, NC))


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
