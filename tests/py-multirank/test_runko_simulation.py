# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import mpi_unittest as mpiut
from mpi4py import MPI
import runko


def single_tile_on_multiple_ranks():
    """
    Here we are explicitly testing that the simulation takes
    correct amount of laps. However, more importantly
    we are implicitly testing that nothing breaks,
    when there less tiles than ranks.
    """

    mpiut.assertEqual(True, 2 <= MPI.COMM_WORLD.Get_size())

    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.n_laps = 7
    config.n_tiles = [1, 1, 1]
    config.n_cells_per_tile = [5, 5, 5]
    config.cfl = 1
    config.field_propagator = "fdtd2"

    tile_grid = runko.TileGrid(config)

    simulation = tile_grid.configure_simulation(config)

    def identity_lap(*_):
        pass

    mpiut.assertEqual(simulation.lap, 0)
    simulation.for_one_lap(identity_lap)
    mpiut.assertEqual(simulation.lap, 1)
    simulation.for_one_lap(identity_lap)
    simulation.for_one_lap(identity_lap)
    mpiut.assertEqual(simulation.lap, 3)
    simulation.for_each_lap(identity_lap)
    mpiut.assertEqual(simulation.lap, 7)


if __name__ == "__main__":
    single_tile_on_multiple_ranks()
