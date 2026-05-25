# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

from mpi4py import MPI


def one_tile_per_rank_cube_size(override_number_of_ranks=None):
    """
    Returns the side length of a cubic tile grid
    such that there is one tile per MPI rank.

    Raises an exception if number of ranks is not a cubic number.

    If override_number_of_ranks (int) is given,
    use it as a number of ranks.
    """

    if override_number_of_ranks:
        ranks = override_number_of_ranks
    else:
        ranks = MPI.COMM_WORLD.size

    candidate = int(round(ranks**(1/3)))
    if candidate**3 != ranks:
        raise RuntimeError(f"Number of ranks {ranks} is not a cubic number!")

    return candidate
