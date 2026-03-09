"""
Multi-rank MPI-IO field snapshot writer tests.

Run with: mpirun -np 4 python test_mpiio_fields.py
"""

import mpi_unittest

import os
import shutil
import tempfile
import numpy as np

from mpi4py import MPI

import runko
from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


HEADER_SIZE = 256
FIELD_NAMES = ["ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz", "rho"]


def read_raw(path):
    """Read a binary MPI-IO field snapshot and return header info + field arrays."""

    with open(path, "rb") as f:
        hdr = f.read(HEADER_SIZE)

    dims = np.frombuffer(hdr[16:28], dtype=np.int32)
    nx, ny, nz = int(dims[0]), int(dims[1]), int(dims[2])

    data = np.fromfile(path, dtype=np.float32, offset=HEADER_SIZE)
    field_size = nx * ny * nz
    fields = {}
    for i, name in enumerate(FIELD_NAMES):
        fields[name] = data[i * field_size : (i + 1) * field_size].reshape(nz, ny, nx)

    return nx, ny, nz, fields


def broadcast_outdir():
    """Create a temp directory on rank 0 and broadcast its path to all ranks."""

    outdir = None
    if rank == 0:
        outdir = tempfile.mkdtemp(prefix="runko-mpiio-test-")
    outdir = comm.bcast(outdir, root=0)
    return outdir


def cleanup_outdir(outdir):
    """Remove the temp directory on rank 0."""

    if rank == 0:
        shutil.rmtree(outdir, ignore_errors=True)


def find_output_file(outdir):
    """Return the path to the first output file in outdir."""

    files = [f for f in os.listdir(outdir) if os.path.isfile(os.path.join(outdir, f))]
    assert len(files) > 0, "No output file was created by MpiioFieldsWriter.write()"
    return os.path.join(outdir, files[0])


# ---------------------------------------------------------------------------
# Test 1: constant fields across all tiles on 4 ranks
# ---------------------------------------------------------------------------
def test_multirank_constant_fields():
    """2x2x2 grid, 8^3 mesh per tile, stride=1, hilbert_curve partitioning.
    Every tile sets constant E=(1,2,3), B=(4,5,6), J=(0,0,0).
    After writing, rank 0 reads the file and verifies all field values."""

    outdir = broadcast_outdir()

    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nt = 1
    config.Nx = 2; config.Ny = 2; config.Nz = 2
    config.NxMesh = 8; config.NyMesh = 8; config.NzMesh = 8
    config.xmin = 0; config.ymin = 0; config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir

    tile_grid = runko.TileGrid(config)

    E_func = lambda x, y, z: (1, 2, 3)
    B_func = lambda x, y, z: (4, 5, 6)
    J_func = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        tile.set_EBJ(E_func, B_func, J_func)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)

    writer = MpiioFieldsWriter(outdir, 2, 8, 2, 8, 2, 8, 1)
    writer.write(tile_grid._corgi_grid, 0)

    # Only rank 0 reads and verifies the file
    if rank == 0:
        path = find_output_file(outdir)
        nx, ny, nz, fields = read_raw(path)

        assert nx == 16, f"Expected nx=16, got {nx}"
        assert ny == 16, f"Expected ny=16, got {ny}"
        assert nz == 16, f"Expected nz=16, got {nz}"

        np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.0, atol=1e-5)
        np.testing.assert_allclose(fields["jx"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["rho"], 0.0, atol=1e-5)

    # Synchronize so all ranks agree the test passed before cleanup
    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 2: tile placement — verify each tile's data lands at the correct
#          global position when written from different ranks
# ---------------------------------------------------------------------------
def test_multirank_tile_placement():
    """2x2x2 grid, 8^3 mesh per tile, stride=1.
    Each tile sets Ex = (i+1) + 10*(j+1) + 100*(k+1) to encode its tile index.
    B, J = (0,0,0).
    After writing, rank 0 reads the global ex array and verifies each 8^3 block
    has the correct encoded value."""

    outdir = broadcast_outdir()

    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nt = 1
    config.Nx = 2; config.Ny = 2; config.Nz = 2
    config.NxMesh = 8; config.NyMesh = 8; config.NzMesh = 8
    config.xmin = 0; config.ymin = 0; config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir

    tile_grid = runko.TileGrid(config)

    B_func = lambda x, y, z: (0, 0, 0)
    J_func = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        ex_val = (idx[0] + 1) + 10 * (idx[1] + 1) + 100 * (idx[2] + 1)
        E_func = lambda x, y, z, v=ex_val: (v, 0, 0)
        tile.set_EBJ(E_func, B_func, J_func)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)

    writer = MpiioFieldsWriter(outdir, 2, 8, 2, 8, 2, 8, 1)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_output_file(outdir)
        nx, ny, nz, fields = read_raw(path)

        assert nx == 16, f"Expected nx=16, got {nx}"
        assert ny == 16, f"Expected ny=16, got {ny}"
        assert nz == 16, f"Expected nz=16, got {nz}"

        ex = fields["ex"]  # shape: (nz, ny, nx) = (16, 16, 16)

        for i_tile in range(2):
            for j_tile in range(2):
                for k_tile in range(2):
                    expected_val = (i_tile + 1) + 10 * (j_tile + 1) + 100 * (k_tile + 1)
                    region = ex[
                        k_tile * 8 : (k_tile + 1) * 8,
                        j_tile * 8 : (j_tile + 1) * 8,
                        i_tile * 8 : (i_tile + 1) * 8,
                    ]
                    np.testing.assert_allclose(
                        region, expected_val, atol=1e-5,
                        err_msg=f"Tile ({i_tile},{j_tile},{k_tile}) Ex mismatch: "
                                f"expected {expected_val}")

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 3: asymmetric mesh dimensions across multiple ranks
# ---------------------------------------------------------------------------
def test_multirank_asymmetric_mesh():
    """2x2x2 grid, 10x11x13 mesh per tile, stride=1.
    Constant fields. Rank 0 verifies global dimensions are (20, 22, 26)
    and field values are correct."""

    outdir = broadcast_outdir()

    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nt = 1
    config.Nx = 2; config.Ny = 2; config.Nz = 2
    config.NxMesh = 10; config.NyMesh = 11; config.NzMesh = 13
    config.xmin = 0; config.ymin = 0; config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir

    tile_grid = runko.TileGrid(config)

    E_func = lambda x, y, z: (7, 8, 9)
    B_func = lambda x, y, z: (10, 11, 12)
    J_func = lambda x, y, z: (0, 0, 0)

    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        tile.set_EBJ(E_func, B_func, J_func)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)

    writer = MpiioFieldsWriter(outdir, 2, 10, 2, 11, 2, 13, 1)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_output_file(outdir)
        nx, ny, nz, fields = read_raw(path)

        assert nx == 20, f"Expected nx=20, got {nx}"
        assert ny == 22, f"Expected ny=22, got {ny}"
        assert nz == 26, f"Expected nz=26, got {nz}"

        for name in FIELD_NAMES:
            assert fields[name].shape == (26, 22, 20), \
                f"Field '{name}' has wrong shape: {fields[name].shape}"

        np.testing.assert_allclose(fields["ex"],  7.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"],  8.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"],  9.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 10.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 11.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 12.0, atol=1e-5)
        np.testing.assert_allclose(fields["jx"],  0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"],  0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"],  0.0, atol=1e-5)
        np.testing.assert_allclose(fields["rho"], 0.0, atol=1e-5)

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


if __name__ == "__main__":
    test_multirank_constant_fields()
    test_multirank_tile_placement()
    test_multirank_asymmetric_mesh()
