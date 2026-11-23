# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import numpy as np
import runko


def make_config(cfl=0.45):
    """electrons at weight 1e-3, positrons, and a massless neutral photon species"""
    config = runko.Configuration(None)
    config.n_tiles = [1, 1, 1]
    config.n_cells_per_tile = [8, 8, 8]
    config.cfl = cfl
    config.field_propagator = "fdtd2"
    config.q0, config.m0, config.w0 = -1, 1, 1e-3
    config.q1, config.m1, config.w1 = 1, 1, 1e-3
    config.q2, config.m2 = 0, 0
    config.particle_pusher = ["boris", "boris", "ballistic"]
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    return config


def make_tile(config=None, idx=(0, 0, 0)):
    config = config or make_config()
    tile = runko.pic.threeD.Tile(idx, config)
    zero = lambda x, y, z: (0, 0, 0)
    tile.set_EBJ(zero, zero, zero)
    return config, tile


class photon_pusher(unittest.TestCase):

    def test_moves_at_c_along_k(self):
        config, tile = make_tile()
        k = [(1e-3, 0, 0), (0, 3.0, 4.0), (-2.0, 2.0, -1.0)]
        w = [0.5, 2.0, 7.0]
        tile.inject(2, [runko.pic.threeD.ParticleState(pos=(4.0, 4.0, 4.0), vel=kk, wgt=ww)
                        for kk, ww in zip(k, w)])
        x0, y0, z0 = tile.get_positions(2)
        tile.push_particles()
        x1, y1, z1 = tile.get_positions(2)
        vx, vy, vz = tile.get_velocities(2)
        for i, kk in enumerate(k):
            n = np.linalg.norm(kk)
            d = np.array([x1[i] - x0[i], y1[i] - y0[i], z1[i] - z0[i]])
            self.assertAlmostEqual(np.linalg.norm(d), config.cfl, places=5)
            np.testing.assert_allclose(d, config.cfl * np.array(kk) / n, rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose([vx[i], vy[i], vz[i]], kk, rtol=1e-6)
        np.testing.assert_allclose(tile.get_weights(2), w, rtol=1e-6)

    def test_lepton_weight_is_the_container_weight(self):
        config, tile = make_tile()
        tile.inject(0, [runko.pic.threeD.ParticleState(pos=(3.0, 3.0, 3.0), vel=(0.1, 0, 0))
                        for _ in range(5)])
        np.testing.assert_allclose(tile.get_weights(0), np.full(5, 1e-3), rtol=1e-6)

    def test_species_without_weight_weighs_one(self):
        config = make_config()
        config.w1 = None
        _, tile = make_tile(config)
        tile.inject(1, [runko.pic.threeD.ParticleState(pos=(3.0, 3.0, 3.0), vel=(0.1, 0, 0))])
        np.testing.assert_allclose(tile.get_weights(1), [1.0])

    def test_photons_deposit_no_current(self):
        """a neutral massless species is skipped by the depositers"""
        for depositer in ("zigzag_1st", "zigzag_1st_atomic"):
            config = make_config()
            config.current_depositer = depositer
            _, tile = make_tile(config)
            rng = np.random.default_rng(2)
            tile.inject(2, [runko.pic.threeD.ParticleState(pos=tuple(rng.uniform(2, 6, 3)),
                                                           vel=tuple(rng.normal(size=3)), wgt=1.0)
                            for _ in range(200)])
            tile.push_particles()
            tile.deposit_current()
            _, _, (jx, jy, jz) = tile.get_EBJ()
            self.assertEqual(float(np.abs(jx).max() + np.abs(jy).max() + np.abs(jz).max()), 0.0)

    def test_neutral_species_needs_the_ballistic_pusher(self):
        config = make_config()
        config.particle_pusher = "boris"
        _, tile = make_tile(config)
        with self.assertRaises(RuntimeError):
            tile.push_particles()

    def test_batch_inject_with_weights(self):
        config, tile = make_tile()
        def pgen(x, y, z):
            n = len(x)
            vel = [np.full(n, 0.01), np.zeros(n), np.zeros(n)]
            return runko.pic.threeD.ParticleStateBatch(pos=[x, y, z], vel=vel, wgt=np.arange(n, dtype=float))
        tile.batch_inject_to_cells(2, pgen)
        w = np.sort(tile.get_weights(2))
        np.testing.assert_allclose(w, np.arange(len(w), dtype=float))
        # leptons ignore a given weight array
        tile.batch_inject_to_cells(0, pgen)
        np.testing.assert_allclose(tile.get_weights(0), np.full(len(w), 1e-3), rtol=1e-6)

    def test_batch_weights_follow_the_stripe_cut(self):
        """entries dropped outside the stripe take their weights with them"""
        config, tile = make_tile()
        def pgen(x, y, z):
            vel = [np.full(len(x), 0.01), np.zeros(len(x)), np.zeros(len(x))]
            return runko.pic.threeD.ParticleStateBatch(pos=[x + 0.5, y, z], vel=vel, wgt=x + 0.5)
        tile.batch_inject_in_x_stripe(2, pgen, 2.0, 4.0)
        xs, _, _ = tile.get_positions(2)
        np.testing.assert_allclose(tile.get_weights(2), xs, rtol=1e-6)
        self.assertTrue(np.all((np.asarray(xs) >= 2.0) & (np.asarray(xs) < 4.0)))

    def test_sort_keeps_weights_paired(self):
        config, tile = make_tile()
        rng = np.random.default_rng(1)
        n = 300
        pos = rng.uniform(0.5, 7.5, (n, 3))
        k = rng.normal(size=(n, 3))
        w = rng.uniform(0.1, 5.0, n)
        tile.inject(2, [runko.pic.threeD.ParticleState(pos=tuple(p), vel=tuple(kk), wgt=ww)
                        for p, kk, ww in zip(pos, k, w)])
        x_before = np.asarray(tile.get_positions(2)[0])
        tile.sort_particles()
        self.assertFalse(np.array_equal(np.asarray(tile.get_positions(2)[0]), x_before))  # reordered
        vx, vy, vz = tile.get_velocities(2)
        ws = tile.get_weights(2)
        self.assertEqual(len(ws), n)
        # weight must follow its own momentum through the permutation
        order = np.argsort(k[:, 0])
        got = np.argsort(vx)
        np.testing.assert_allclose(np.asarray(ws)[got], w[order], rtol=1e-6)

    def test_weights_survive_tile_communication(self):
        """photons crossing into the neighbouring tile keep their weights"""
        import tempfile
        config = make_config()
        config.tile_partitioning = "catepillar_track"
        config.catepillar_track_length = 1
        config.n_laps = 1
        config.n_tiles = [2, 1, 1]
        config.io_outdir = tempfile.mkdtemp(prefix="runko-photon-comm-")
        tile_grid = runko.TileGrid(config)
        rng = np.random.default_rng(3)
        n = 200
        w = rng.uniform(0.1, 5.0, n)
        for idx in tile_grid.local_tile_indices():
            _, tile = make_tile(config, idx)
            if tuple(idx) == (0, 0, 0):
                # at the +x edge of the first tile, flying +x: one push moves them into the second
                tile.inject(2, [runko.pic.threeD.ParticleState(pos=(7.8, 2.0 + 4.0 * rng.random(), 2.0 + 4.0 * rng.random()),
                                                              vel=(1.0, 0.0, 0.0), wgt=ww) for ww in w])
            tile_grid.add_tile(tile, idx)
        simulation = tile_grid.configure_simulation(config)
        def lap(x):
            x.prtcl_push()
            x.comm_external(runko.tools.comm_mode.pic_particle)
            x.comm_local(runko.tools.comm_mode.pic_particle)
        simulation.for_each_lap(lap)
        got = {}
        for tile in simulation.local_tiles():
            xs = np.asarray(tile.get_positions(2)[0]); ws = np.asarray(tile.get_weights(2))
            got[tuple(tile.index)] = (xs, ws)
        x1, w1 = got[(1, 0, 0)]
        self.assertEqual(len(w1), n)
        self.assertEqual(len(got[(0, 0, 0)][1]), 0)
        self.assertTrue(np.all(x1 >= 8.0))
        np.testing.assert_allclose(np.sort(w1), np.sort(w), rtol=1e-6)


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
