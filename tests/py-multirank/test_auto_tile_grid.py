# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

# Assumed to be run with 8 ranks.

import mpi_unittest
import runko


def deduce_one_tile_per_rank_cube():
    N = runko.one_tile_per_rank_cube_size()
    mpi_unittest.assertEqual(2, N)


if __name__ == "__main__":
    deduce_one_tile_per_rank_cube()
