import unittest
import tempfile
import shutil

import numpy as np

import runko
from tests.mpiio_test_helpers import (
    FIELD_NAMES, field_names, make_config, make_pic_config,
    setup_grid_with_fields, setup_pic_grid_with_particles,
    find_output_file, write_and_read,
)
from runko.mpiio_reader import MAGIC, read_header


class TestMpiioFieldsWriterConstruction(unittest.TestCase):

    # TODO re-use the setup, teardown defined below and combine this test
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
        self.assertEqual(hdr["num_fields"], 11)

        for i, expected_name in enumerate(FIELD_NAMES):
            self.assertEqual(hdr["field_names"][i], expected_name)

    def test_single_tile_field_values(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Constant E=(1,2,3), B=(4,5,6), J=(7,8,9).
        Read back binary data and verify all field arrays have correct values.
        Runs for both write() and write_collective()."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (7, 8, 9))

        for collective in (False, True):
            with self.subTest(collective=collective):
                _, fields = write_and_read(
                    tile_grid, self.outdir, config, collective=collective)

                np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
                np.testing.assert_allclose(fields["ey"], 2.0, atol=1e-5)
                np.testing.assert_allclose(fields["ez"], 3.0, atol=1e-5)
                np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)
                np.testing.assert_allclose(fields["by"], 5.0, atol=1e-5)
                np.testing.assert_allclose(fields["bz"], 6.0, atol=1e-5)
                np.testing.assert_allclose(fields["jx"], 7.0, atol=1e-5)
                np.testing.assert_allclose(fields["jy"], 8.0, atol=1e-5)
                np.testing.assert_allclose(fields["jz"], 9.0, atol=1e-5)
                np.testing.assert_allclose(fields["n0"], 0.0, atol=1e-5)
                np.testing.assert_allclose(fields["n1"], 0.0, atol=1e-5)

    def test_asymmetric_mesh(self):
        """1x1x1 grid, 10x11x13 mesh, stride=1. All dims different.
        Set constant fields and verify output shape is (13, 11, 10)."""

        config = make_config(NxMesh=10, NyMesh=11, NzMesh=13, outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (7, 8, 9),
            lambda x, y, z: (10, 11, 12),
            lambda x, y, z: (13, 14, 15))

        hdr, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["nx"], 10)
        self.assertEqual(hdr["ny"], 11)
        self.assertEqual(hdr["nz"], 13)

        for name in FIELD_NAMES:
            self.assertEqual(fields[name].shape, (13, 11, 10),
                             f"Field '{name}' has wrong shape: {fields[name].shape}")

        np.testing.assert_allclose(fields["ex"],  7.0, atol=1e-5)
        np.testing.assert_allclose(fields["ey"],  8.0, atol=1e-5)
        np.testing.assert_allclose(fields["ez"],  9.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 10.0, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 11.0, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 12.0, atol=1e-5)
        np.testing.assert_allclose(fields["jx"], 13.0, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], 14.0, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], 15.0, atol=1e-5)

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

    def test_stride_equals_mesh_size(self):
        """stride=NxMesh produces a single output cell per tile.
        E/B are subsampled at cell [0,0,0].
        J is volume-summed over all NxMesh^3 cells."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (7, 8, 9))

        stride = 8  # equals NxMesh
        hdr, fields = write_and_read(
            tile_grid, self.outdir, config, stride=stride)

        self.assertEqual(hdr["nx"], 1)
        self.assertEqual(hdr["ny"], 1)
        self.assertEqual(hdr["nz"], 1)

        # E/B: subsampled value at cell [0,0,0]
        np.testing.assert_allclose(fields["ex"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.0, atol=1e-5)

        # J: volume-summed over 8^3 = 512 cells
        volume = stride ** 3
        np.testing.assert_allclose(fields["jx"], volume * 7.0, atol=1e-3)
        np.testing.assert_allclose(fields["jy"], volume * 8.0, atol=1e-3)
        np.testing.assert_allclose(fields["jz"], volume * 9.0, atol=1e-3)

    def test_multi_tile_placement(self):
        """2x2x1 grid, 8^3 mesh per tile, stride=1, catepillar_track_length=1.
        Set different constant Ex in each tile. Verify global array has values
        at correct positions. Global output is 16x16x8.
        Runs for both write() and write_collective()."""

        config = make_config(Nx=2, Ny=2, Nz=1, outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.emf.threeD.Tile(idx, config)

            base = (idx[0] + 1) + 10 * (idx[1] + 1)

            E_func = lambda x, y, z, v=base: (v,       v + 100, v + 200)
            B_func = lambda x, y, z, v=base: (v + 300, v + 400, v + 500)
            J_func = lambda x, y, z, v=base: (v + 600, v + 700, v + 800)
            tile.set_EBJ(E_func, B_func, J_func)
            tile_grid.add_tile(tile, idx)

        _ = tile_grid.configure_simulation(config)

        field_offsets = {
            "ex": 0,   "ey": 100, "ez": 200,
            "bx": 300, "by": 400, "bz": 500,
            "jx": 600, "jy": 700, "jz": 800,
        }

        for collective in (False, True):
            with self.subTest(collective=collective):
                hdr, fields = write_and_read(
                    tile_grid, self.outdir, config, collective=collective)

                self.assertEqual(hdr["nx"], 16)
                self.assertEqual(hdr["ny"], 16)
                self.assertEqual(hdr["nz"], 8)

                for fname, offset in field_offsets.items():
                    for i_tile in range(2):
                        for j_tile in range(2):
                            expected_val = (i_tile + 1) + 10 * (j_tile + 1) + offset
                            region = fields[fname][
                                :,
                                j_tile * 8:(j_tile + 1) * 8,
                                i_tile * 8:(i_tile + 1) * 8]
                            np.testing.assert_allclose(
                                region, expected_val, atol=1e-5,
                                err_msg=f"Tile ({i_tile},{j_tile},0) {fname} mismatch")

    def test_reader_roundtrip(self):
        """Write with known values, read with pytools reader, compare."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1.5, 2.5, 3.5),
            lambda x, y, z: (4.5, 5.5, 6.5),
            lambda x, y, z: (7.5, 8.5, 9.5))

        _, fields = write_and_read(tile_grid, self.outdir, config)

        self.assertIsInstance(fields, dict)
        self.assertEqual(len(fields), 11)

        for name in FIELD_NAMES:
            self.assertIn(name, fields)
            self.assertEqual(fields[name].shape, (8, 8, 8))

        np.testing.assert_allclose(fields["ex"], 1.5, atol=1e-5)
        np.testing.assert_allclose(fields["ey"], 2.5, atol=1e-5)
        np.testing.assert_allclose(fields["ez"], 3.5, atol=1e-5)
        np.testing.assert_allclose(fields["bx"], 4.5, atol=1e-5)
        np.testing.assert_allclose(fields["by"], 5.5, atol=1e-5)
        np.testing.assert_allclose(fields["bz"], 6.5, atol=1e-5)
        np.testing.assert_allclose(fields["jx"], 7.5, atol=1e-5)
        np.testing.assert_allclose(fields["jy"], 8.5, atol=1e-5)
        np.testing.assert_allclose(fields["jz"], 9.5, atol=1e-5)
        np.testing.assert_allclose(fields["n0"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["n1"], 0.0, atol=1e-5)

    def test_spatially_varying_fields(self):
        """1x1x1 grid, 8^3 mesh, stride=1. Set E=(x, y, z) where x,y,z are
        cell-center coordinates (0.5, 1.5, ..., 7.5). Verify per-cell data
        ordering. This catches x/y/z transposition bugs that constant-field
        tests miss. Runs for both write() and write_collective()."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (x,    y,    z),
            lambda x, y, z: (x+10, y+10, z+10),
            lambda x, y, z: (x+20, y+20, z+20))

        for collective in (False, True):
            with self.subTest(collective=collective):
                _, fields = write_and_read(
                    tile_grid, self.outdir, config, collective=collective)

                # Output is C-order (nz, ny, nx): fields[k, j, i]
                # Yee staggering in set_EBJ:
                #   Ex(x+0.5,y,z)     → ex[k,j,i] = i + 0.5
                #   Ey(x,y+0.5,z)     → ey[k,j,i] = j + 0.5
                #   Ez(x,y,z+0.5)     → ez[k,j,i] = k + 0.5
                #   Bx(x,y+0.5,z+0.5) → bx[k,j,i] = i + 10   (x is integer)
                #   By(x+0.5,y,z+0.5) → by[k,j,i] = j + 10   (y is integer)
                #   Bz(x+0.5,y+0.5,z) → bz[k,j,i] = k + 10   (z is integer)
                #   Jx(x+0.5,y,z)     → jx[k,j,i] = i + 20.5 (same stagger as E)
                #   Jy(x,y+0.5,z)     → jy[k,j,i] = j + 20.5
                #   Jz(x,y,z+0.5)     → jz[k,j,i] = k + 20.5
                for k in range(8):
                    for j in range(8):
                        for i in range(8):
                            self.assertAlmostEqual(
                                fields["ex"][k, j, i], i + 0.5, places=4,
                                msg=f"ex[{k},{j},{i}] != {i+0.5}")
                            self.assertAlmostEqual(
                                fields["ey"][k, j, i], j + 0.5, places=4,
                                msg=f"ey[{k},{j},{i}] != {j+0.5}")
                            self.assertAlmostEqual(
                                fields["ez"][k, j, i], k + 0.5, places=4,
                                msg=f"ez[{k},{j},{i}] != {k+0.5}")
                            self.assertAlmostEqual(
                                fields["bx"][k, j, i], i + 10, places=4,
                                msg=f"bx[{k},{j},{i}] != {i+10}")
                            self.assertAlmostEqual(
                                fields["by"][k, j, i], j + 10, places=4,
                                msg=f"by[{k},{j},{i}] != {j+10}")
                            self.assertAlmostEqual(
                                fields["bz"][k, j, i], k + 10, places=4,
                                msg=f"bz[{k},{j},{i}] != {k+10}")
                            self.assertAlmostEqual(
                                fields["jx"][k, j, i], i + 20.5, places=4,
                                msg=f"jx[{k},{j},{i}] != {i+20.5}")
                            self.assertAlmostEqual(
                                fields["jy"][k, j, i], j + 20.5, places=4,
                                msg=f"jy[{k},{j},{i}] != {j+20.5}")
                            self.assertAlmostEqual(
                                fields["jz"][k, j, i], k + 20.5, places=4,
                                msg=f"jz[{k},{j},{i}] != {k+20.5}")

    def test_header_completeness(self):
        """Verify all header fields are correct with asymmetric grid/mesh/stride."""

        config = make_config(Nx=2, Ny=3, Nz=1,
                             NxMesh=10, NyMesh=11, NzMesh=13,
                             outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (7, 8, 9))

        stride = 2
        hdr, _ = write_and_read(tile_grid, self.outdir, config, stride=stride)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["version"], 3)
        self.assertEqual(hdr["header_size"], 512)
        self.assertEqual(hdr["dtype_size"], 4)
        self.assertEqual(hdr["stride"], stride)
        self.assertEqual(hdr["lap"], 0)
        self.assertEqual(hdr["Nx"], 2)
        self.assertEqual(hdr["Ny"], 3)
        self.assertEqual(hdr["Nz"], 1)
        self.assertEqual(hdr["NxMesh"], 10)
        self.assertEqual(hdr["NyMesh"], 11)
        self.assertEqual(hdr["NzMesh"], 13)
        # Output dims: Nx*(NxMesh//stride) = 2*5=10, 3*5=15, 1*6=6
        self.assertEqual(hdr["nx"], 10)
        self.assertEqual(hdr["ny"], 15)
        self.assertEqual(hdr["nz"], 6)

    def test_nonzero_lap(self):
        """Write with lap=5, verify header and filename contain lap number."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (0, 0, 0))

        import os
        hdr, _ = write_and_read(tile_grid, self.outdir, config, lap=5)

        self.assertEqual(hdr["lap"], 5)

        # Verify filename contains the lap number
        files = os.listdir(self.outdir)
        self.assertTrue(any("5" in f for f in files),
                        f"No file with lap=5 in {files}")


class TestMpiioFieldsDensity(unittest.TestCase):

    # TODO combine this class with the one above

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_uniform_density_ppc1(self):
        """PIC tile, 1 particle per cell per species, stride=1.
        n0 and n1 should both be 1.0 everywhere.
        Runs for both write() and write_collective()."""

        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        for collective in (False, True):
            with self.subTest(collective=collective):
                _, fields = write_and_read(
                    tile_grid, self.outdir, config, collective=collective)
                np.testing.assert_allclose(fields["n0"], 1.0, atol=1e-5)
                np.testing.assert_allclose(fields["n1"], 1.0, atol=1e-5)

    def test_uniform_density_ppc4(self):
        """PIC tile, 4 particles per cell per species, stride=1.
        n0 and n1 should both be 4.0 everywhere."""

        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=4)

        _, fields = write_and_read(tile_grid, self.outdir, config)
        np.testing.assert_allclose(fields["n0"], 4.0, atol=1e-5)
        np.testing.assert_allclose(fields["n1"], 4.0, atol=1e-5)

    def test_density_with_stride(self):
        """PIC tile, ppc=1, stride=2.
        Each output cell sums stride^3=8 full-res cells.
        With 1 particle per cell, n0 and n1 should be 8.0."""

        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        stride = 2
        _, fields = write_and_read(tile_grid, self.outdir, config, stride=stride)
        volume = stride ** 3
        np.testing.assert_allclose(fields["n0"], volume * 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["n1"], volume * 1.0, atol=1e-5)

    def test_emf_tile_density_zero(self):
        """EMF-only tile (no particles). n0 and n1 should be 0."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (0, 0, 0))

        _, fields = write_and_read(tile_grid, self.outdir, config)
        np.testing.assert_allclose(fields["n0"], 0.0, atol=1e-5)
        np.testing.assert_allclose(fields["n1"], 0.0, atol=1e-5)


class TestMpiioFieldsMultiSpecies(unittest.TestCase):

    # combine this class with the first 

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_emf_only_zero_species(self):
        """EMF-only tile with nspecies=0. Only 9 EMF fields, no density fields."""

        config = make_config(outdir=self.outdir)
        tile_grid = setup_grid_with_fields(
            config,
            lambda x, y, z: (1, 2, 3),
            lambda x, y, z: (4, 5, 6),
            lambda x, y, z: (0, 0, 0))

        hdr, fields = write_and_read(tile_grid, self.outdir, config, nspecies=0)
        self.assertEqual(hdr["num_fields"], 9)
        self.assertEqual(len(fields), 9)
        self.assertNotIn("n0", fields)

    def test_three_species(self):
        """PIC tile with 3 species. n0, n1, n2 should all be ppc."""

        nsp = 3
        config = make_pic_config(outdir=self.outdir, nspecies=nsp)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(
            config, zero, zero, zero, ppc=1, nspecies=nsp)

        hdr, fields = write_and_read(
            tile_grid, self.outdir, config, nspecies=nsp)

        self.assertEqual(hdr["num_fields"], 9 + nsp)
        self.assertEqual(hdr["field_names"], field_names(nsp))
        for i in range(nsp):
            np.testing.assert_allclose(fields[f"n{i}"], 1.0, atol=1e-5)

    def test_five_species(self):
        """PIC tile with 5 species (max). n0-n4 should all be ppc."""

        nsp = 5
        config = make_pic_config(outdir=self.outdir, nspecies=nsp)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(
            config, zero, zero, zero, ppc=2, nspecies=nsp)

        hdr, fields = write_and_read(
            tile_grid, self.outdir, config, nspecies=nsp)

        self.assertEqual(hdr["num_fields"], 14)
        for i in range(nsp):
            np.testing.assert_allclose(fields[f"n{i}"], 2.0, atol=1e-5)

    def test_one_species(self):
        """PIC tile with 1 species. Only n0 field written."""

        nsp = 1
        config = make_pic_config(outdir=self.outdir, nspecies=nsp)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(
            config, zero, zero, zero, ppc=1, nspecies=nsp)

        hdr, fields = write_and_read(
            tile_grid, self.outdir, config, nspecies=nsp)

        self.assertEqual(hdr["num_fields"], 10)
        np.testing.assert_allclose(fields["n0"], 1.0, atol=1e-5)
        self.assertNotIn("n1", fields)


if __name__ == "__main__":
    unittest.main()
