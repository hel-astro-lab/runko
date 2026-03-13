import unittest
import tempfile
import shutil
import os

import numpy as np

import runko
from runko_cpp_bindings.emf.threeD import MpiioParticlesWriter
from runko.mpiio_prtcl_reader import MAGIC, read_prtcl_header, read_prtcl_snapshot
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

    # TODO combine this class with the others; re-use setUp and tearDown

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

    def test_position_dependent_velocity(self):
        """Inject particles with vel = (x+0.5, y+0.5, z+0.5), verify roundtrip."""
        P = runko.pic.threeD.ParticleState
        config = make_pic_config(outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            def gen(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5),
                          vel=(x + 0.5, y + 0.5, z + 0.5))]
            for sp in range(2):
                tile.inject_to_each_cell(sp, gen)
            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        total = 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        # Each particle's velocity should match its position
        np.testing.assert_allclose(fields["ux"], fields["x"], atol=1e-4)
        np.testing.assert_allclose(fields["uy"], fields["y"], atol=1e-4)
        np.testing.assert_allclose(fields["uz"], fields["z"], atol=1e-4)

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

    def test_linear_field_interpolation(self):
        """Linearly varying E=(x,y,z). Interior particles see correct
        interpolated fields."""
        config = make_pic_config(outdir=self.outdir)

        E_func = lambda x, y, z: (x, y, z)
        B_func = lambda x, y, z: (0, 0, 0)
        J_func = lambda x, y, z: (0, 0, 0)

        tile_grid = setup_pic_grid_with_particles(
            config, E_func, B_func, J_func, ppc=1)

        total = 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        x, y, z = fields["x"], fields["y"], fields["z"]
        interior = ((x > 1.5) & (x < 6.5) &
                    (y > 1.5) & (y < 6.5) &
                    (z > 1.5) & (z < 6.5))

        self.assertGreater(np.sum(interior), 0, "No interior particles found")

        # Yee staggering: Ex at (x+0.5, y, z). For a linear field E_func=(x,y,z),
        # Ex = x+0.5 at the staggered position. The interpolator should recover
        # the field at the particle position. Start with generous tolerance.
        np.testing.assert_allclose(
            fields["ex"][interior], x[interior], atol=1.0)
        np.testing.assert_allclose(
            fields["ey"][interior], y[interior], atol=1.0)
        np.testing.assert_allclose(
            fields["ez"][interior], z[interior], atol=1.0)

    def test_two_species_separate_files(self):
        """Write species 0 and 1 separately, verify different files and
        that both species have correct velocities."""
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

        h0 = read_prtcl_header(find_prtcl_file(self.outdir, 0))
        h1 = read_prtcl_header(find_prtcl_file(self.outdir, 1))
        f0 = read_prtcl_snapshot(find_prtcl_file(self.outdir, 0))
        f1 = read_prtcl_snapshot(find_prtcl_file(self.outdir, 1))

        self.assertEqual(h0["species"], 0)
        self.assertEqual(h1["species"], 1)

        np.testing.assert_allclose(f0["ux"], 1.0, atol=1e-5)
        np.testing.assert_allclose(f0["uy"], 0.0, atol=1e-5)
        np.testing.assert_allclose(f1["ux"], 0.0, atol=1e-5)
        np.testing.assert_allclose(f1["uy"], 2.0, atol=1e-5)

    def test_multi_tile_positions(self):
        """2x2x1 grid, 8^3 mesh, ppc=1. Encode tile index in velocity,
        then verify each tile's particles occupy a distinct NxMesh-sized
        spatial block and that global positions are non-overlapping."""
        config = make_pic_config(Nx=2, Ny=2, Nz=1, outdir=self.outdir)
        P = runko.pic.threeD.ParticleState

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Encode tile index as velocity code
            code = float(100 * idx[0] + 10 * idx[1] + idx[2])
            def gen(x, y, z, _c=code):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(_c, 0, 0))]
            for species in range(2):
                tile.inject_to_each_cell(species, gen)
            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        total = 4 * 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        self.assertEqual(hdr["n_prtcls"], total)

        # Group particles by velocity code and verify spatial properties
        tile_codes = sorted(set(fields["ux"].astype(int)))
        self.assertEqual(len(tile_codes), 4, "Expected 4 distinct tiles")

        tile_x_ranges = []
        for code_val in tile_codes:
            mask = fields["ux"] == float(code_val)
            n_in_tile = np.sum(mask)
            self.assertEqual(n_in_tile, 8**3,
                             f"Tile code={code_val}: expected {8**3} particles, "
                             f"got {n_in_tile}")

            x_min, x_max = fields["x"][mask].min(), fields["x"][mask].max()
            y_min, y_max = fields["y"][mask].min(), fields["y"][mask].max()

            # Each tile should span exactly NxMesh cells (range = 7.0)
            self.assertAlmostEqual(x_max - x_min, 7.0, places=2,
                                   msg=f"Tile code={code_val} x range wrong")
            self.assertAlmostEqual(y_max - y_min, 7.0, places=2,
                                   msg=f"Tile code={code_val} y range wrong")

            tile_x_ranges.append((x_min, x_max))

        # Verify tiles at different grid x-indices have non-overlapping x positions
        x_starts = sorted(r[0] for r in tile_x_ranges)
        self.assertGreater(x_starts[-1] - x_starts[0], 8 - 1,
                           "Tiles do not occupy distinct x positions")

    def test_multi_tile_velocities(self):
        """2x2x1 grid with different velocities per tile.
        Verifies tile_buf_offset correctly interleaves particles."""
        config = make_pic_config(Nx=2, Ny=2, Nz=1, outdir=self.outdir)
        P = runko.pic.threeD.ParticleState

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Encode tile index in velocity
            vx = float(idx[0] + 1)
            vy = float(idx[1] + 1)

            def gen(x, y, z, _vx=vx, _vy=vy):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(_vx, _vy, 0))]
            for species in range(2):
                tile.inject_to_each_cell(species, gen)

            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        total = 4 * 8**3
        hdr, fields = prtcl_write_and_read(
            tile_grid, self.outdir, n_prtcls=total * 10)

        # Each tile should contribute particles with its encoded velocity
        for i_tile in range(2):
            for j_tile in range(2):
                x_lo, x_hi = i_tile * 8, (i_tile + 1) * 8
                y_lo, y_hi = j_tile * 8, (j_tile + 1) * 8
                in_tile = ((fields["x"] >= x_lo) & (fields["x"] < x_hi) &
                           (fields["y"] >= y_lo) & (fields["y"] < y_hi))
                expected_vx = float(i_tile + 1)
                expected_vy = float(j_tile + 1)
                np.testing.assert_allclose(
                    fields["ux"][in_tile], expected_vx, atol=1e-5,
                    err_msg=f"Tile ({i_tile},{j_tile}) ux mismatch")
                np.testing.assert_allclose(
                    fields["uy"][in_tile], expected_vy, atol=1e-5,
                    err_msg=f"Tile ({i_tile},{j_tile}) uy mismatch")

    def test_empty_species(self):
        """Species 0 has particles, species 1 has zero.
        Writing species 1 should produce a valid file with n_prtcls=0."""
        config = make_pic_config(outdir=self.outdir, nspecies=2)
        P = runko.pic.threeD.ParticleState

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Only inject species 0
            def gen(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(1, 0, 0))]
            tile.inject_to_each_cell(0, gen)
            # Species 1 gets no particles

            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        writer = MpiioParticlesWriter(self.outdir, 1000, 1)
        writer.write(tile_grid._corgi_grid, 0)

        path = find_prtcl_file(self.outdir, 1)
        hdr = read_prtcl_header(path)
        self.assertEqual(hdr["n_prtcls"], 0)
        self.assertEqual(hdr["species"], 1)

    def test_header_completeness(self):
        """Verify all header fields are correct."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        hdr, _ = prtcl_write_and_read(tile_grid, self.outdir, n_prtcls=10000)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["version"], 3)
        self.assertEqual(hdr["header_size"], 512)
        self.assertEqual(hdr["dtype_size"], 4)
        self.assertEqual(hdr["species"], 0)
        self.assertEqual(hdr["lap"], 0)
        self.assertEqual(hdr["num_fields"], 12)

    def test_nonzero_lap(self):
        """Write with lap=3, verify header and filename."""
        config = make_pic_config(outdir=self.outdir)
        zero = lambda x, y, z: (0, 0, 0)
        tile_grid = setup_pic_grid_with_particles(config, zero, zero, zero, ppc=1)

        hdr, _ = prtcl_write_and_read(tile_grid, self.outdir, n_prtcls=10000, lap=3)

        self.assertEqual(hdr["lap"], 3)
        files = os.listdir(self.outdir)
        self.assertTrue(any("3" in f for f in files),
                        f"No file with lap=3 in {files}")


if __name__ == "__main__":
    unittest.main()
