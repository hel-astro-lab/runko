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
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    return config


def make_tile(config=None):
    config = config or make_config()
    tile = runko.pic.threeD.Tile((0, 0, 0), config)
    zero = lambda x, y, z: (0, 0, 0)
    tile.set_EBJ(zero, zero, zero)
    return config, tile


class photon_pusher(unittest.TestCase):

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


if __name__ == "__main__":
    _ = runko.actions.RuntimeInstance()
    unittest.main()
