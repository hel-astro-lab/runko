import mpi_unittest
import numpy as np

import runko

def create_test_grid():
    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.Nt = 3
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


def emf_halos_are_initialized():
    """
    TileGrid with only emf tiles, which have been initialized
    to have some constant E and B vector fields (curl(E or B) = 0).
    If halos have not ben initialized, then curl(E or B) != 0
    and pushing fields will not keep them constant.
    """

    conf, tile_grid = create_test_grid()

    E0 = lambda x, y, z: (1, 2, 3)
    B0 = lambda x, y, z: (4, 5, 6)
    J0 = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.Tile(idx, conf)
        tile.set_EBJ(E0, B0, J0)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    def assertConstantFields():
        asserts = []
        A = lambda x: mpi_unittest.assertEqualDeferred(True, x)

        for tile in simulation.local_tiles():
            (Ex, Ey, Ez), (Bx, By, Bz), _ = tile.get_EBJ()

            asserts.append(A(np.all(Ex == 1)))
            asserts.append(A(np.all(Ey == 2)))
            asserts.append(A(np.all(Ez == 3)))
            asserts.append(A(np.all(Bx == 4)))
            asserts.append(A(np.all(By == 5)))
            asserts.append(A(np.all(Bz == 6)))

        mpi_unittest.assertDeferredResults(asserts)
        mpi_unittest.assertEqual(True, len(asserts) >= 1)

    # Initial fields are constant.
    assertConstantFields()

    def f(local_tile, communicate, *_):
        EBmodes = (runko.comm_mode.emf_E, runko.comm_mode.emf_B)

        communicate.virtual_tile_sync(*EBmodes)
        communicate.pairwise_moore(*EBmodes)

        local_tile.push_half_b()
        local_tile.push_e()

    # Fields stay constant even after pushing them.
    simulation.for_one_lap(f)
    assertConstantFields()

    simulation.for_each_lap(f)
    assertConstantFields()


if __name__ == "__main__":
    virtual_emf_tiles()
    emf_halos_are_initialized()
