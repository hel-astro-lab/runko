import unittest
import tempfile
import shutil
import os

import numpy as np

import runko


HEADER_SIZE = 256
MAGIC = 0x524E4B4F
FIELD_NAMES = ["ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz", "rho"]


def read_raw(path):
    """Read a binary MPI-IO field snapshot and return header info + field arrays."""

    with open(path, "rb") as f:
        hdr = f.read(HEADER_SIZE)

    magic = np.frombuffer(hdr[0:4], dtype=np.uint32)[0]
    num_fields = int(np.frombuffer(hdr[12:16], dtype=np.uint32)[0])
    dims = np.frombuffer(hdr[16:28], dtype=np.int32)
    nx, ny, nz = int(dims[0]), int(dims[1]), int(dims[2])

    # Field names: 10 names, each 16 chars null-padded, at bytes 64-224
    names = []
    for i in range(num_fields):
        raw = hdr[64 + i * 16 : 64 + (i + 1) * 16]
        names.append(raw.split(b"\x00")[0].decode("ascii"))

    data = np.fromfile(path, dtype=np.float32, offset=HEADER_SIZE)
    field_size = nx * ny * nz
    fields = {}
    for i, name in enumerate(FIELD_NAMES):
        fields[name] = data[i * field_size : (i + 1) * field_size].reshape(nz, ny, nx)

    return magic, nx, ny, nz, num_fields, names, fields


def make_config(Nx=1, Ny=1, Nz=1, NxMesh=8, NyMesh=8, NzMesh=8, outdir=None):
    """Create a standard test configuration."""

    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.Nx = Nx
    config.Ny = Ny
    config.Nz = Nz
    config.NxMesh = NxMesh
    config.NyMesh = NyMesh
    config.NzMesh = NzMesh
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.Nt = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir
    return config


def setup_grid_with_fields(config, E_func, B_func, J_func):
    """Create a TileGrid, populate tiles with fields, and configure simulation."""

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        tile.set_EBJ(E_func, B_func, J_func)
        tile_grid.add_tile(tile, idx)

    simulation = tile_grid.configure_simulation(config)
    return tile_grid, simulation


class TestMpiioFieldsWriterConstruction(unittest.TestCase):

    def test_construction(self):
        """Create MpiioFieldsWriter and verify no crash."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        outdir = tempfile.mkdtemp()
        try:
            writer = MpiioFieldsWriter(outdir, 1, 8, 1, 8, 1, 8, 1)
        finally:
            shutil.rmtree(outdir)


class TestMpiioFieldsSingleTile(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_single_tile_write(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Constant E=(1,2,3), B=(4,5,6), J=(0,0,0).
        Verify file exists and header is correct."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        config = make_config(outdir=self.outdir)
        E_func = lambda x, y, z: (1, 2, 3)
        B_func = lambda x, y, z: (4, 5, 6)
        J_func = lambda x, y, z: (0, 0, 0)
        tile_grid, simulation = setup_grid_with_fields(config, E_func, B_func, J_func)

        writer = MpiioFieldsWriter(self.outdir, 1, 8, 1, 8, 1, 8, 1)
        writer.write(tile_grid._corgi_grid, 0)

        # Find the output file
        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        self.assertTrue(len(output_files) > 0,
                        "No output file was created by MpiioFieldsWriter.write()")

        # Read and verify header
        output_path = os.path.join(self.outdir, output_files[0])
        magic, nx, ny, nz, num_fields, names, fields = read_raw(output_path)

        self.assertEqual(magic, MAGIC)
        self.assertEqual(nx, 8)
        self.assertEqual(ny, 8)
        self.assertEqual(nz, 8)
        self.assertEqual(num_fields, 10)

        # Verify field names in header
        for i, expected_name in enumerate(FIELD_NAMES):
            self.assertEqual(names[i], expected_name)


    def test_single_tile_field_values(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Constant E=(1,2,3), B=(4,5,6), J=(0,0,0).
        Read back binary data and verify all field arrays have correct values."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        config = make_config(outdir=self.outdir)
        E_func = lambda x, y, z: (1, 2, 3)
        B_func = lambda x, y, z: (4, 5, 6)
        J_func = lambda x, y, z: (0, 0, 0)
        tile_grid, simulation = setup_grid_with_fields(config, E_func, B_func, J_func)

        writer = MpiioFieldsWriter(self.outdir, 1, 8, 1, 8, 1, 8, 1)
        writer.write(tile_grid._corgi_grid, 0)

        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        output_path = os.path.join(self.outdir, output_files[0])
        magic, nx, ny, nz, num_fields, names, fields = read_raw(output_path)

        # E fields: constant values 1, 2, 3
        np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.0, atol=1e-5)

        # B fields: constant values 4, 5, 6
        np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.0, atol=1e-5)

        # J fields: constant zero
        np.testing.assert_allclose(fields["jx"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], 0.0, atol=1e-5)

        # rho: always zero (not stored in EBJ)
        np.testing.assert_allclose(fields["rho"], 0.0, atol=1e-5)


    def test_asymmetric_mesh(self):
        """1x1x1 grid, 10x11x13 mesh, stride=1. All dims different.
        Set constant fields and verify output shape is (13, 11, 10)."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        config = make_config(NxMesh=10, NyMesh=11, NzMesh=13, outdir=self.outdir)
        E_func = lambda x, y, z: (7, 8, 9)
        B_func = lambda x, y, z: (10, 11, 12)
        J_func = lambda x, y, z: (0, 0, 0)
        tile_grid, simulation = setup_grid_with_fields(config, E_func, B_func, J_func)

        writer = MpiioFieldsWriter(self.outdir, 1, 10, 1, 11, 1, 13, 1)
        writer.write(tile_grid._corgi_grid, 0)

        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        output_path = os.path.join(self.outdir, output_files[0])
        magic, nx, ny, nz, num_fields, names, fields = read_raw(output_path)

        # Verify dimensions
        self.assertEqual(nx, 10)
        self.assertEqual(ny, 11)
        self.assertEqual(nz, 13)

        # Verify shape: (nz, ny, nx)
        for name in FIELD_NAMES:
            self.assertEqual(fields[name].shape, (13, 11, 10),
                             f"Field '{name}' has wrong shape: {fields[name].shape}")

        # Verify constant field values
        np.testing.assert_allclose(fields["ex"], 7.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 8.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 9.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 10.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 11.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 12.0, atol=1e-5)


    def test_stride(self):
        """1x1x1 grid, 8^3 mesh, stride=2. Set constant fields.
        Verify output shape is 4x4x4. For E/B verify subsampled values.
        For J verify volume-averaged values."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        config = make_config(outdir=self.outdir)
        E_func = lambda x, y, z: (1, 2, 3)
        B_func = lambda x, y, z: (4, 5, 6)
        J_func = lambda x, y, z: (1, 1, 1)
        tile_grid, simulation = setup_grid_with_fields(config, E_func, B_func, J_func)

        stride = 2
        writer = MpiioFieldsWriter(self.outdir, 1, 8, 1, 8, 1, 8, stride)
        writer.write(tile_grid._corgi_grid, 0)

        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        output_path = os.path.join(self.outdir, output_files[0])
        magic, nx, ny, nz, num_fields, names, fields = read_raw(output_path)

        # Verify strided dimensions
        self.assertEqual(nx, 4)
        self.assertEqual(ny, 4)
        self.assertEqual(nz, 4)

        for name in FIELD_NAMES:
            self.assertEqual(fields[name].shape, (4, 4, 4),
                             f"Field '{name}' has wrong shape: {fields[name].shape}")

        # E/B: subsampled from constant field, values should be preserved
        np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.0, atol=1e-5)

        # J: volume-averaged from constant field.
        # Each strided cell averages stride^3 = 8 identical values of 1.0.
        # The existing FieldsWriter sums (does not divide), so the expected
        # value is stride^3 * J = 8.0 for each component.
        volume = stride ** 3
        np.testing.assert_allclose(fields["jx"], volume * 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], volume * 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], volume * 1.0, atol=1e-5)


class TestMpiioFieldsMultiTile(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_multi_tile_placement(self):
        """2x2x1 grid, 8^3 mesh per tile, stride=1, catepillar_track_length=1.
        Set different constant Ex in each tile. Verify global array has values
        at correct positions. Global output is 16x16x8."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        config = make_config(Nx=2, Ny=2, Nz=1, outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.emf.threeD.Tile(idx, config)
            # Set Ex to a value that encodes tile position: (i+1) + 10*(j+1)
            ex_val = (idx[0] + 1) + 10 * (idx[1] + 1)
            E_func = lambda x, y, z, v=ex_val: (v, 0, 0)
            B_func = lambda x, y, z: (0, 0, 0)
            J_func = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(E_func, B_func, J_func)
            tile_grid.add_tile(tile, idx)

        simulation = tile_grid.configure_simulation(config)

        writer = MpiioFieldsWriter(self.outdir, 2, 8, 2, 8, 1, 8, 1)
        writer.write(tile_grid._corgi_grid, 0)

        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        output_path = os.path.join(self.outdir, output_files[0])
        magic, nx, ny, nz, num_fields, names, fields = read_raw(output_path)

        # Verify global dimensions
        self.assertEqual(nx, 16)
        self.assertEqual(ny, 16)
        self.assertEqual(nz, 8)

        ex = fields["ex"]  # shape: (8, 16, 16)

        # Check that each tile region has the correct constant value.
        # Tile (i_tile, j_tile, 0) occupies global x in [i_tile*8 : (i_tile+1)*8],
        #                                      global y in [j_tile*8 : (j_tile+1)*8].
        # In the (nz, ny, nx) array: ex[:, j_tile*8:(j_tile+1)*8, i_tile*8:(i_tile+1)*8]
        for i_tile in range(2):
            for j_tile in range(2):
                expected_val = (i_tile + 1) + 10 * (j_tile + 1)
                region = ex[:, j_tile * 8:(j_tile + 1) * 8,
                               i_tile * 8:(i_tile + 1) * 8]
                np.testing.assert_allclose(
                    region, expected_val, atol=1e-5,
                    err_msg=f"Tile ({i_tile},{j_tile},0) Ex mismatch")


class TestMpiioFieldsReaderRoundtrip(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_reader_roundtrip(self):
        """Write with known values, read with pytools reader, compare."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter
        from pytools.mpiio_reader import read_field_snapshot

        config = make_config(outdir=self.outdir)
        E_func = lambda x, y, z: (1.5, 2.5, 3.5)
        B_func = lambda x, y, z: (4.5, 5.5, 6.5)
        J_func = lambda x, y, z: (0, 0, 0)
        tile_grid, simulation = setup_grid_with_fields(config, E_func, B_func, J_func)

        writer = MpiioFieldsWriter(self.outdir, 1, 8, 1, 8, 1, 8, 1)
        writer.write(tile_grid._corgi_grid, 0)

        output_files = [f for f in os.listdir(self.outdir)
                        if os.path.isfile(os.path.join(self.outdir, f))]
        output_path = os.path.join(self.outdir, output_files[0])

        fields = read_field_snapshot(output_path)

        self.assertIsInstance(fields, dict)
        self.assertEqual(len(fields), 10)

        # Verify shape
        for name in FIELD_NAMES:
            self.assertIn(name, fields)
            self.assertEqual(fields[name].shape, (8, 8, 8))

        # Verify values
        np.testing.assert_allclose(fields["ex"], 1.5, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.5, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.5, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.5, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.5, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.5, atol=1e-5)
        np.testing.assert_allclose(fields["jx"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["rho"], 0.0, atol=1e-5)


if __name__ == "__main__":
    unittest.main()
