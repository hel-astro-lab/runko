import mpi_unittest
import runko

def create_test_grid():
    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.Nt = 1
    config.Nx = 2
    config.Ny = 2
    config.Nz = 2
    config.NxMesh = 10
    config.NyMesh = 11
    config.NzMesh = 13
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"

    return config, runko.TileGrid(config)


def virtual_emf_tiles():
    """
    Here we test that after grid.configure_simulation(...)
    each virtual tile has bee initialized to emf tile,
    as they are next to them.
    """

    conf, tile_grid = create_test_grid()

    mpi_unittest.assertEqual(tile_grid.initialized_from_restart_file(), False)

    local_idx = list(tile_grid.local_tile_indices())
    mpi_unittest.assertNotEqual(len(local_idx), 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.Tile(idx, conf)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    vtiles = list(simulation.virtual_tiles())
    mpi_unittest.assertNotEqual(len(vtiles), 0)

    asserts = []
    for vtile in vtiles:
        asserts.append(mpi_unittest.assertEqualDeferred(type(vtile), runko.emf.Tile))

    mpi_unittest.assertDeferredResults(asserts)


if __name__ == "__main__":
    virtual_emf_tiles()
