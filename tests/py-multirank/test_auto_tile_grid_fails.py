# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

# Assumed to be run with 3 ranks.

import mpi_unittest
import runko


def deduce_one_tile_per_rank_cube():
    ok = False
    try:
        N = runko.one_tile_per_rank_cube_size()
    except:
        ok = True
    mpi_unittest.assertEqual(ok, True)


if __name__ == "__main__":
    deduce_one_tile_per_rank_cube()
