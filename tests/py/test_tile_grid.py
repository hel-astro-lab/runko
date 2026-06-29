# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import runko

class tile_grid(unittest.TestCase):
    def setUp(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "hilbert_curve"
        conf.n_tiles = [2, 2, 2]
        conf.n_cells_per_tile = [10, 11, 13]

        self.test_grid = runko.TileGrid(conf)

    def test_mandatory_configuration_variables(self):
        empty_conf = runko.Configuration(None)

        with self.assertRaises(RuntimeError):
            runko.TileGrid(empty_conf)

    def test_nonsense_tile_partitioning(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "nonsense"
        conf.n_tiles = [2, 2, 2]
        conf.n_cells_per_tile = [10, 11, 13]

        with self.assertRaisesRegex(RuntimeError, r"tile_partitioning"):
            runko.TileGrid(conf)

    def test_catepillar_track_requires_length(self):
        conf = runko.Configuration(None)
        conf.tile_partitioning = "catepillar_track"
        conf.n_tiles = [2, 2, 2]
        conf.n_cells_per_tile = [10, 11, 13]

        with self.assertRaisesRegex(RuntimeError, r"catepillar_track_length"):
            runko.TileGrid(conf)

if __name__ == "__main__":
    unittest.main()
