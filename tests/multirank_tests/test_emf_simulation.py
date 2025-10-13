import mpi_unittest
import numpy as np
import tempfile

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
    config.outdir = tempfile.mkdtemp(prefix="runko-emf-test-output")

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
        tile = runko.emf.threeD.Tile(idx, conf)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    vtiles = list(simulation.virtual_tiles())
    mpi_unittest.assertNotEqual(len(vtiles), 0)

    asserts = []
    for vtile in vtiles:
        asserts.append(mpi_unittest.assertEqualDeferred(type(vtile), runko.emf.threeD.Tile))

    mpi_unittest.assertDeferredResults(asserts)


def emf_communication():
    """
    TileGrid with only emf tiles, which have been initialized
    to have some constant E and B vector fields (curl(E or B) = 0).
    Initially the halo regions of each tile have not been initialized,
    which means that curl(E or B) != 0 next to halo
    and pushing fields will not keep them constant.

    In order to initialize the halos one has to do virtual tile sync
    for each tile and then pariwise moore communication.
    """

    conf, tile_grid = create_test_grid()

    E0 = lambda x, y, z: (1, 2, 3)
    B0 = lambda x, y, z: (4, 5, 6)
    J0 = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, conf)
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
        EBmodes = (runko.tools.comm_mode.emf_E, runko.tools.comm_mode.emf_B)

        communicate.virtual_tile_sync(*EBmodes)
        communicate.pairwise_moore(*EBmodes)

        local_tile.push_half_b()
        local_tile.push_e()

    # Fields stay constant even after pushing them.
    simulation.for_one_lap(f)
    assertConstantFields()

    simulation.for_each_lap(f)
    assertConstantFields()


def emf_J_exchange():
    """
    TileGrid with only emf tiles, which have been initialized
    to have some constant J0. Then after virtual tile sync
    normal J pairwaise moore communication and another virtual tile sync,
    the halo regions also have the current J0 (*).

    Now when pairwise moore J exchange is executed,
    all non-halo J should be either J0, 2 * J0, 4 * J0 or 8 * J0.
    It is fiddly to check exactly that right cells have correct amounts,
    se we only check that there are correct amount of each.

    (*) The first virtual tile sync sets J0 to non-halo regions of virtual tiles.
    Pairwise moore sets J0 to halo regions of local tiles.
    The second virtual tile sync sets J0 to non-halo regions of virtual tiles.
    """

    conf, tile_grid = create_test_grid()

    E = lambda x, y, z: (0, 0, 0)
    B = lambda x, y, z: (0, 0, 0)

    J0 = np.array((1, 2, 3))

    def J(x, y, z):
        return J0

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, conf)
        tile.set_EBJ(E, B, J)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(conf)

    def f(tile, communicate, *_):
        communicate.virtual_tile_sync(runko.tools.comm_mode.emf_J)
        communicate.pairwise_moore(runko.tools.comm_mode.emf_J)
        communicate.virtual_tile_sync(runko.tools.comm_mode.emf_J)
        communicate.pairwise_moore(runko.tools.comm_mode.emf_J_exchange)

    simulation.for_one_lap(f)


    nx, ny, nz = conf.NxMesh, conf.NyMesh, conf.NzMesh
    halo_width = 3
    twohalo = 2 * halo_width

    expected_num_of_J0 = (nx - twohalo) * (ny - twohalo) * (nz - twohalo)

    exp_num_of_2J0_a = 2 * (nx - twohalo) * (ny - twohalo)
    exp_num_of_2J0_b = 2 * (ny - twohalo) * (nz - twohalo)
    exp_num_of_2J0_c = 2 * (nx - twohalo) * (nz - twohalo)
    expected_num_of_2J0 = exp_num_of_2J0_a + exp_num_of_2J0_b + exp_num_of_2J0_c

    exp_num_of_4J0_a = 4 * (nx - twohalo) * halo_width
    exp_num_of_4J0_b = 4 * (ny - twohalo) * halo_width
    exp_num_of_4J0_c = 4 * (nz - twohalo) * halo_width
    expected_num_of_4J0 = exp_num_of_4J0_a + exp_num_of_4J0_b + exp_num_of_4J0_c

    expected_num_of_8J0 = 8 * halo_width * halo_width

    asserts = []

    for tile in simulation.local_tiles():
        _, _, (Jx, Jy, Jz) = tile.get_EBJ()

        # Assume that np.unique returns sorted unique values.
        unique_Jx, counts_Jx = np.unique(Jx, return_counts=True)
        unique_Jy, counts_Jy = np.unique(Jy, return_counts=True)
        unique_Jz, counts_Jz = np.unique(Jz, return_counts=True)

        mpi_unittest.assertEqualDeferred(len(unique_Jx), 4)
        mpi_unittest.assertEqualDeferred(len(unique_Jy), 4)
        mpi_unittest.assertEqualDeferred(len(unique_Jz), 4)

        mpi_unittest.assertEqualDeferred(unique_Jx[0], J0[0])
        mpi_unittest.assertEqualDeferred(unique_Jx[1], 2 * J0[0])
        mpi_unittest.assertEqualDeferred(unique_Jx[2], 4 * J0[0])
        mpi_unittest.assertEqualDeferred(unique_Jx[3], 8 * J0[0])

        mpi_unittest.assertEqualDeferred(unique_Jy[0], J0[1])
        mpi_unittest.assertEqualDeferred(unique_Jy[1], 2 * J0[1])
        mpi_unittest.assertEqualDeferred(unique_Jy[2], 4 * J0[1])
        mpi_unittest.assertEqualDeferred(unique_Jy[3], 8 * J0[1])

        mpi_unittest.assertEqualDeferred(unique_Jz[0], J0[2])
        mpi_unittest.assertEqualDeferred(unique_Jz[1], 2 * J0[2])
        mpi_unittest.assertEqualDeferred(unique_Jz[2], 4 * J0[2])
        mpi_unittest.assertEqualDeferred(unique_Jz[3], 8 * J0[2])

        mpi_unittest.assertEqualDeferred(counts_Jx[0], expected_num_of_J0)
        mpi_unittest.assertEqualDeferred(counts_Jx[1], expected_num_of_2J0)
        mpi_unittest.assertEqualDeferred(counts_Jx[2], expected_num_of_4J0)
        mpi_unittest.assertEqualDeferred(counts_Jx[3], expected_num_of_8J0)

        mpi_unittest.assertEqualDeferred(counts_Jy[0], expected_num_of_J0)
        mpi_unittest.assertEqualDeferred(counts_Jy[1], expected_num_of_2J0)
        mpi_unittest.assertEqualDeferred(counts_Jy[2], expected_num_of_4J0)
        mpi_unittest.assertEqualDeferred(counts_Jy[3], expected_num_of_8J0)

        mpi_unittest.assertEqualDeferred(counts_Jz[0], expected_num_of_J0)
        mpi_unittest.assertEqualDeferred(counts_Jz[1], expected_num_of_2J0)
        mpi_unittest.assertEqualDeferred(counts_Jz[2], expected_num_of_4J0)
        mpi_unittest.assertEqualDeferred(counts_Jz[3], expected_num_of_8J0)


    mpi_unittest.assertDeferredResults(asserts)


if __name__ == "__main__":
    virtual_emf_tiles()
    emf_communication()
    emf_J_exchange()
