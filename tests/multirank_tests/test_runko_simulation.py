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
    config.Nt = 7
    config.Nx = 1
    config.Ny = 1
    config.Nz = 1
    config.NxMesh = 5
    config.NyMesh = 5
    config.NzMesh = 5
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"

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
