import unittest
import tempfile
import shutil

import numpy as np

import runko
from tests.mpiio_test_helpers import (
    FIELD_NAMES, make_config, setup_grid_with_fields,
    find_output_file, write_and_read,
)
from pytools.mpiio_reader import MAGIC, read_header


class TestMpiioFieldsWriterConstruction(unittest.TestCase):

    def test_construction(self):
        """Create MpiioFieldsWriter and verify no crash."""

        from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter

        outdir = tempfile.mkdtemp()
        try:
            writer = MpiioFieldsWriter(outdir, 1, 8, 1, 8, 1, 8, 1)
        finally:
            shutil.rmtree(outdir)


class TestMpiioFieldsWriter(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_single_tile_write(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Constant E=(1,2,3), B=(4,5,6), J=(0,0,0).
        Verify file exists and header is correct."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (0, 0, 0))

        hdr, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["nx"], 8)
        self.assertEqual(hdr["ny"], 8)
        self.assertEqual(hdr["nz"], 8)
        self.assertEqual(hdr["num_fields"], 10)

        for i, expected_name in enumerate(FIELD_NAMES):
            self.assertEqual(hdr["field_names"][i], expected_name)

    def test_single_tile_field_values(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Constant E=(1,2,3), B=(4,5,6), J=(0,0,0).
        Read back binary data and verify all field arrays have correct values."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (0, 0, 0))

        _, fields = write_and_read(tile_grid, self.outdir, config)

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

    def test_asymmetric_mesh(self):
        """1x1x1 grid, 10x11x13 mesh, stride=1. All dims different.
        Set constant fields and verify output shape is (13, 11, 10)."""

        config = make_config(NxMesh=10, NyMesh=11, NzMesh=13, outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (7, 8, 9),
            lambda x, y, z: (10, 11, 12),
            lambda x, y, z: (0, 0, 0))

        hdr, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["nx"], 10)
        self.assertEqual(hdr["ny"], 11)
        self.assertEqual(hdr["nz"], 13)

        for name in FIELD_NAMES:
            self.assertEqual(fields[name].shape, (13, 11, 10),
                             f"Field '{name}' has wrong shape: {fields[name].shape}")

        np.testing.assert_allclose(fields["ex"], 7.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 8.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 9.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 10.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 11.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 12.0, atol=1e-5)

    def test_stride(self):
        """1x1x1 grid, 8^3 mesh, stride=2. Set constant fields.
        Verify output shape is 4x4x4. For E/B verify subsampled values.
        For J verify volume-summed values."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (1, 1, 1))

        stride = 2
        hdr, fields = write_and_read(tile_grid, self.outdir, config, stride=stride)

        self.assertEqual(hdr["nx"], 4)
        self.assertEqual(hdr["ny"], 4)
        self.assertEqual(hdr["nz"], 4)

        for name in FIELD_NAMES:
            self.assertEqual(fields[name].shape, (4, 4, 4),
                             f"Field '{name}' has wrong shape: {fields[name].shape}")

        np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.0, atol=1e-5)

        # J: volume-sum from constant field.
        # Each strided cell sums stride^3 = 8 identical values of 1.0.
        volume = stride ** 3
        np.testing.assert_allclose(fields["jx"], volume * 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], volume * 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], volume * 1.0, atol=1e-5)

    def test_multi_tile_placement(self):
        """2x2x1 grid, 8^3 mesh per tile, stride=1, catepillar_track_length=1.
        Set different constant Ex in each tile. Verify global array has values
        at correct positions. Global output is 16x16x8."""

        config = make_config(Nx=2, Ny=2, Nz=1, outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.emf.threeD.Tile(idx, config)
            ex_val = (idx[0] + 1) + 10 * (idx[1] + 1)
            E_func = lambda x, y, z, v=ex_val: (v, 0, 0)
            B_func = lambda x, y, z: (0, 0, 0)
            J_func = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(E_func, B_func, J_func)
            tile_grid.add_tile(tile, idx)

        _ = tile_grid.configure_simulation(config)

        hdr, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["nx"], 16)
        self.assertEqual(hdr["ny"], 16)
        self.assertEqual(hdr["nz"], 8)

        ex = fields["ex"]  # shape: (8, 16, 16)

        for i_tile in range(2):
            for j_tile in range(2):
                expected_val = (i_tile + 1) + 10 * (j_tile + 1)
                region = ex[:, j_tile * 8:(j_tile + 1) * 8,
                               i_tile * 8:(i_tile + 1) * 8]
                np.testing.assert_allclose(
                    region, expected_val, atol=1e-5,
                    err_msg=f"Tile ({i_tile},{j_tile},0) Ex mismatch")

    def test_reader_roundtrip(self):
        """Write with known values, read with pytools reader, compare."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1.5, 2.5, 3.5),
            lambda x, y, z: (4.5, 5.5, 6.5),
            lambda x, y, z: (0, 0, 0))

        _, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertIsInstance(fields, dict)
        self.assertEqual(len(fields), 10)

        for name in FIELD_NAMES:
            self.assertIn(name, fields)
            self.assertEqual(fields[name].shape, (8, 8, 8))

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
