# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import numpy as np
import runko


def make_config(cfl=0.45):
    """electrons, positrons, and a massless neutral photon species"""
    config = runko.Configuration(None)
    config.n_tiles = [1, 1, 1]
    config.n_cells_per_tile = [8, 8, 8]
    config.cfl = cfl
    config.field_propagator = "fdtd2"
    config.q0, config.m0 = -1, 1
    config.q1, config.m1 = 1, 1
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
        tile.inject(2, [runko.pic.threeD.ParticleState(pos=(4.0, 4.0, 4.0), vel=kk) for kk in k])
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

    def test_photons_deposit_no_current(self):
        """a neutral massless species is skipped by the depositers"""
        for depositer in ("zigzag_1st", "zigzag_1st_atomic"):
            config = make_config()
            config.current_depositer = depositer
            _, tile = make_tile(config)
            rng = np.random.default_rng(2)
            tile.inject(2, [runko.pic.threeD.ParticleState(pos=tuple(rng.uniform(2, 6, 3)),
                                                           vel=tuple(rng.normal(size=3)))
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


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
