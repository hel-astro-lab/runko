# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

from mpi4py import MPI


def one_tile_per_rank_cube_size():
    """
    Returns the side length of a cubic tile grid
    such that there is one tile per MPI rank.

    Raises an exception if number of ranks is not a cubic number.
    """
    ranks = MPI.COMM_WORLD.size

    if int(ranks**(1/3))**3 != ranks:
        raise RuntimeError("Number of ranks is not a cubic number!")

    return int(ranks**(1/3))
