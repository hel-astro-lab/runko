# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import runko


class auto_tile_grid(unittest.TestCase):
    def test_one_tile_per_rank_cube_size(self):
        ok_ranks = [1, 8, 27, 64, 125, 512]
        expected_sizes = [1, 2, 3, 4, 5, 8]
        for rank, size in zip(ok_ranks, expected_sizes):
            n = runko.one_tile_per_rank_cube_size(rank)
            self.assertEqual(n, size)

        not_ok_ranks = [2, 7, 28, 63, 65, 126, 511, 513]
        for rank in not_ok_ranks:
            with self.assertRaises(Exception):
                runko.one_tile_per_rank_cube_size(rank)


if __name__ == "__main__":
    unittest.main()
