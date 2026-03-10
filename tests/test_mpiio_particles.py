import unittest
import tempfile
import shutil
import os

import numpy as np

import runko
from runko_cpp_bindings.emf.threeD import MpiioParticlesWriter
from pytools.mpiio_prtcl_reader import MAGIC, read_prtcl_header, read_prtcl_snapshot
from tests.mpiio_test_helpers import (
    make_pic_config, setup_pic_grid_with_particles, PRTCL_FIELD_NAMES,
)


def find_prtcl_file(outdir, species=0):
    """Return path to the particle snapshot file for given species."""
    files = [f for f in os.listdir(outdir)
             if f.startswith(f"prtcls_{species}_") and f.endswith(".bin")]
    assert len(files) > 0, f"No particle output file for species {species}"
    return os.path.join(outdir, files[0])


def prtcl_write_and_read(tile_grid, outdir, n_prtcls, species=0, lap=0):
    """Write a particle snapshot and read it back.

    Returns (header_dict, fields_dict).
    """
    writer = MpiioParticlesWriter(outdir, n_prtcls, species)
    writer.write(tile_grid._corgi_grid, lap)
    path = find_prtcl_file(outdir, species)
    hdr = read_prtcl_header(path)
    fields = read_prtcl_snapshot(path)
    return hdr, fields


class TestMpiioParticlesWriterConstruction(unittest.TestCase):

    def test_construction(self):
        """Create MpiioParticlesWriter and verify no crash."""
        outdir = tempfile.mkdtemp()
        try:
            writer = MpiioParticlesWriter(outdir, 100, 0)
        finally:
            shutil.rmtree(outdir)


class TestMpiioParticlesWriter(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_header_fields(self):
        """Verify header has correct magic, version, field names."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        hdr, _ = prtcl_write_and_read(tile_grid, self.outdir, n_prtcls=100)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["num_fields"], 12)
        self.assertEqual(hdr["field_names"], PRTCL_FIELD_NAMES)
        self.assertEqual(hdr["species"], 0)
        self.assertGreater(hdr["n_prtcls"], 0)

    def test_all_particles_when_n_exceeds_total(self):
        """When n_prtcls > total particles, all particles should be written."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        total_particles = 8 * 8 * 8  # ppc=1, 8^3 mesh
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total_particles * 10)

        self.assertEqual(hdr["n_prtcls"], total_particles)
        for name in PRTCL_FIELD_NAMES:
            self.assertEqual(len(fields[name]), total_particles)

    def test_sampling_reduces_count(self):
        """When n_prtcls < total, output should have fewer particles."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=2)

        total_particles = 2 * 8 * 8 * 8  # ppc=2, 8^3 mesh = 1024
        n_target = 100
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=n_target)

        self.assertGreater(hdr["n_prtcls"], 0)
        self.assertLess(hdr["n_prtcls"], total_particles)
        # Allow some tolerance due to integer rounding
        self.assertAlmostEqual(hdr["n_prtcls"], n_target, delta=n_target * 0.5)

    def test_position_values_global(self):
        """Verify sampled positions are in global coordinates.
        1x1x1 grid with mins at (0,0,0), so positions should be in [0, NxMesh)."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        total = 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        # All positions should be within [0, NxMesh=8)
        for coord in ("x", "y", "z"):
            self.assertTrue(np.all(fields[coord] >= 0),
                            f"{coord} has negative values")
            self.assertTrue(np.all(fields[coord] < 8),
                            f"{coord} has values >= NxMesh")

    def test_velocity_values(self):
        """Inject particles with known velocities and verify output."""
        P = runko.pic.threeD.ParticleState
        config = make_pic_config(outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Inject particles with known velocity
            def gen(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(1.0, 2.0, 3.0))]
            tile.inject_to_each_cell(0, gen)

            # Species 1 with different velocity
            def gen2(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(4.0, 5.0, 6.0))]
            tile.inject_to_each_cell(1, gen2)

            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        total = 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10, species=0)

        np.testing.assert_allclose(fields["ux"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["uy"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["uz"], 3.0, atol=1e-5)

    def test_uniform_field_interpolation(self):
        """Set uniform E=(1,2,3), B=(4,5,6). Verify interpolated fields in output."""
        config = make_pic_config(outdir=self.outdir)

        E_func = lambda x, y, z: (1, 2, 3)
        B_func = lambda x, y, z: (4, 5, 6)
        J_func = lambda x, y, z: (0, 0, 0)

        tile_grid = setup_pic_grid_with_particles(config, E_func, B_func, J_func, ppc=1)

        total = 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        # Interior particles (away from halo boundaries) should see uniform fields.
        # Near edges the staggered Yee interpolation may sample halo cells.
        # Filter to particles at interior positions (1 < pos < 7).
        x, y, z = fields["x"], fields["y"], fields["z"]
        interior = (x > 1) & (x < 7) & (y > 1) & (y < 7) & (z > 1) & (z < 7)

        self.assertGreater(np.sum(interior), 0, "No interior particles found")

        np.testing.assert_allclose(fields["ex"][interior], 1.0, atol=1e-4)
        np.testing.assert_allclose(fields["ey"][interior], 2.0, atol=1e-4)
        np.testing.assert_allclose(fields["ez"][interior], 3.0, atol=1e-4)
        np.testing.assert_allclose(fields["bx"][interior], 4.0, atol=1e-4)
        np.testing.assert_allclose(fields["by"][interior], 5.0, atol=1e-4)
        np.testing.assert_allclose(fields["bz"][interior], 6.0, atol=1e-4)

    def test_two_species_separate_files(self):
        """Write species 0 and 1 separately, verify different files."""
        config = make_pic_config(outdir=self.outdir)
        P = runko.pic.threeD.ParticleState

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Species 0: velocity (1, 0, 0)
            def gen0(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(1.0, 0, 0))]
            tile.inject_to_each_cell(0, gen0)

            # Species 1: velocity (0, 2, 0)
            def gen1(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(0, 2.0, 0))]
            tile.inject_to_each_cell(1, gen1)

            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        total = 8**3

        # Write species 0
        w0 = MpiioParticlesWriter(self.outdir, total * 10, 0)
        w0.write(tile_grid._corgi_grid, 0)

        # Write species 1
        w1 = MpiioParticlesWriter(self.outdir, total * 10, 1)
        w1.write(tile_grid._corgi_grid, 0)

        f0 = read_prtcl_snapshot(find_prtcl_file(self.outdir, 0))
        f1 = read_prtcl_snapshot(find_prtcl_file(self.outdir, 1))

        np.testing.assert_allclose(f0["ux"], 1.0, atol=1e-5)
        np.testing.assert_allclose(f0["uy"], 0.0, atol=1e-5)
        np.testing.assert_allclose(f1["ux"], 0.0, atol=1e-5)
        np.testing.assert_allclose(f1["uy"], 2.0, atol=1e-5)


if __name__ == "__main__":
    unittest.main()
