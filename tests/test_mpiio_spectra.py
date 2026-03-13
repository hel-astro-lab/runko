import unittest
import tempfile
import shutil
import os

import numpy as np

import runko
from runko_cpp_bindings.emf.threeD import MpiioSpectraWriter
from runko.mpiio_spectra_reader import (
    read_spectra_header, read_spectra_snapshot,
    u_bin_edges, u_bin_centers, beta_bin_edges, beta_bin_centers,
)
from runko.mpiio_constants import MAGIC
from tests.mpiio_test_helpers import make_pic_config


def find_spectra_file(outdir):
    """Return path to the first spectra output file in outdir."""
    files = [f for f in os.listdir(outdir)
             if f.startswith("pspectra_") and f.endswith(".bin")]
    assert len(files) > 0, "No spectra output file"
    return os.path.join(outdir, files[0])


def make_grid_with_velocity(config, vx=0.0, vy=0.0, vz=0.0, ppc=1, nspecies=2):
    """Create a PIC grid where every cell has ppc particles with given velocity."""
    P = runko.pic.threeD.ParticleState

    def particle_gen(x, y, z):
        return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(vx, vy, vz))
                for _ in range(ppc)]

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, config)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)
        for species in range(nspecies):
            tile.inject_to_each_cell(species, particle_gen)
        tile_grid.add_tile(tile, idx)

    _ = tile_grid.configure_simulation(config)
    return tile_grid


def spectra_write_and_read(tile_grid, outdir, config, nbins=50,
                           umin=1e-2, umax=1e2, stride=1, nspecies=2, lap=0):
    """Write spectra and read back via Python reader."""
    writer = MpiioSpectraWriter(
        outdir,
        config.Nx, config.NxMesh,
        config.Ny, config.NyMesh,
        config.Nz, config.NzMesh,
        stride, nbins, umin, umax, nspecies)
    writer.write(tile_grid._corgi_grid, lap)
    path = find_spectra_file(outdir)
    hdr = read_spectra_header(path)
    fields = read_spectra_snapshot(path)
    return hdr, fields


class TestMpiioSpectraWriterConstruction(unittest.TestCase):

    def test_construction(self):
        """Create MpiioSpectraWriter and verify no crash."""
        outdir = tempfile.mkdtemp()
        try:
            writer = MpiioSpectraWriter(outdir, 1, 8, 1, 8, 1, 8, 1, 50, 1e-2, 1e2)
        finally:
            shutil.rmtree(outdir)


class TestMpiioSpectraWriter(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def test_header_fields(self):
        """Verify header magic, version, field names, nbins, umin, umax."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=2.0, vz=3.0)
        hdr, _ = spectra_write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["num_fields"], 8)  # 2 species * 4
        self.assertEqual(hdr["nbins"], 50)
        self.assertAlmostEqual(hdr["umin"], 1e-2, places=5)
        self.assertAlmostEqual(hdr["umax"], 1e2, places=3)
        self.assertEqual(hdr["field_names"][:4], ["s0_u", "s0_bx", "s0_by", "s0_bz"])
        self.assertEqual(hdr["field_names"][4:8], ["s1_u", "s1_bx", "s1_by", "s1_bz"])

    def test_output_shape(self):
        """Verify output shape is (Nz, Ny, nx_global, nbins)."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=2.0, vz=3.0)
        hdr, fields = spectra_write_and_read(tile_grid, self.outdir, config, nbins=30)

        for name, arr in fields.items():
            self.assertEqual(arr.shape, (1, 1, 8, 30))  # Nz=1, Ny=1, nx=8, nbins=30

    def test_stride_reduces_x(self):
        """stride=2 should produce nx_global = Nx * (NxMesh/2)."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=1.0, vz=1.0)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, stride=2)

        self.assertEqual(hdr["nx"], 4)  # 8/2 = 4
        for name, arr in fields.items():
            self.assertEqual(arr.shape[2], 4)

    def test_known_velocity_u_spectrum(self):
        """Inject particles with ux=3, uy=4, uz=0 => |u|=5.
        All particles should land in the bin containing u=5."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=3.0, vy=4.0, vz=0.0)
        nbins = 100
        umin, umax = 1e-2, 1e4
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=nbins, umin=umin, umax=umax)

        edges = u_bin_edges(hdr)
        expected_bin = np.searchsorted(edges, 5.0) - 1

        s0_u = fields["s0_u"]
        total = s0_u.sum(axis=(0, 1, 2))  # sum over nz, ny, nx -> (nbins,)

        total_particles = 8**3  # ppc=1, 8^3 mesh
        self.assertEqual(total[expected_bin], total_particles)

        # All other bins should be zero
        mask = np.ones(nbins, dtype=bool)
        mask[expected_bin] = False
        np.testing.assert_array_equal(total[mask], 0)

    def test_known_velocity_beta_spectrum(self):
        """Inject particles with ux=1, uy=0, uz=0.
        gamma = sqrt(1+1) = sqrt(2), beta_x = 1/sqrt(2) ~ 0.707.
        Check that the peak lands in the correct linear bin."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=0.0, vz=0.0)
        nbins = 100
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=nbins, umin=1e-2, umax=1e2)

        edges = beta_bin_edges(hdr)
        beta_x = 1.0 / np.sqrt(2.0)
        expected_bin = np.searchsorted(edges, beta_x) - 1

        s0_bx = fields["s0_bx"]
        total = s0_bx.sum(axis=(0, 1, 2))

        total_particles = 8**3
        self.assertEqual(total[expected_bin], total_particles)

        # NOTE: beta_y and beta_z are tested in test_beta_y_and_z below

    def test_beta_particle_count_conserved(self):
        """Sum over all beta_x bins should equal total particles
        (all particles have well-defined beta)."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=2.0, vz=3.0)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=100, umin=1e-4, umax=1e4)

        total_particles = 8**3
        # beta_x is always in [-1, 1], so all particles are counted
        bx_total = fields["s0_bx"].sum()
        self.assertAlmostEqual(bx_total, total_particles, delta=1)

    def test_u_out_of_range_clamped_to_boundary_bins(self):
        """Particles outside [umin, umax] are clamped: above to last bin,
        below to first bin. Total count is always conserved."""
        config = make_pic_config(outdir=self.outdir)
        nbins = 50

        for label, vx, umin, umax, expected_bin in [
            ("above umax", 100.0, 1e-2, 10.0, nbins - 1),
            ("below umin", 0.001, 1.0, 1e4,   0),
        ]:
            with self.subTest(case=label):
                tile_grid = make_grid_with_velocity(config, vx=vx, vy=0.0, vz=0.0)
                hdr, fields = spectra_write_and_read(
                    tile_grid, self.outdir, config,
                    nbins=nbins, umin=umin, umax=umax)

                total_particles = 8**3
                total = fields["s0_u"].sum(axis=(0, 1, 2))

                self.assertEqual(total[expected_bin], total_particles)
                self.assertAlmostEqual(total.sum(), total_particles, delta=1)

    def test_u_particle_count_conserved(self):
        """Sum over all u bins should equal total particles
        (out-of-range particles go to last bin)."""
        config = make_pic_config(outdir=self.outdir)
        # ux=1, uy=2, uz=3 => |u| = sqrt(14) ~ 3.74
        tile_grid = make_grid_with_velocity(config, vx=1.0, vy=2.0, vz=3.0)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=100, umin=1e-4, umax=1e4)

        total_particles = 8**3
        u_total = fields["s0_u"].sum()
        self.assertAlmostEqual(u_total, total_particles, delta=1)

    def test_bin_edge_helpers(self):
        """Verify bin edge helper functions."""
        hdr = {"nbins": 10, "umin": 0.01, "umax": 100.0}

        u_edges = u_bin_edges(hdr)
        self.assertEqual(len(u_edges), 11)
        self.assertAlmostEqual(u_edges[0], 0.01)
        self.assertAlmostEqual(u_edges[-1], 100.0)

        u_cents = u_bin_centers(hdr)
        self.assertEqual(len(u_cents), 10)

        b_edges = beta_bin_edges(hdr)
        self.assertEqual(len(b_edges), 11)
        self.assertAlmostEqual(b_edges[0], -1.0)
        self.assertAlmostEqual(b_edges[-1], 1.0)

        b_cents = beta_bin_centers(hdr)
        self.assertEqual(len(b_cents), 10)
        self.assertAlmostEqual(b_cents[5], 0.1, places=5)  # midpoint ~ 0

    def test_multi_tile_x_resolved(self):
        """2x1x1 grid, each tile with different velocity.
        Verify the x-resolved spectra have peaks at different bins for each tile."""
        config = make_pic_config(Nx=2, Ny=1, Nz=1, outdir=self.outdir)
        P = runko.pic.threeD.ParticleState

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            # Tile 0: |u|=1, tile 1: |u|=50
            vx = 1.0 if idx[0] == 0 else 50.0

            def gen(x, y, z, _vx=vx):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(_vx, 0, 0))]
            for species in range(2):
                tile.inject_to_each_cell(species, gen)
            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        nbins = 100
        umin, umax = 1e-2, 1e4
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=nbins, umin=umin, umax=umax)

        # Output shape: (Nz=1, Ny=1, nx=16, nbins)
        self.assertEqual(fields["s0_u"].shape, (1, 1, 16, nbins))

        edges = u_bin_edges(hdr)
        bin_u1 = np.searchsorted(edges, 1.0) - 1
        bin_u50 = np.searchsorted(edges, 50.0) - 1
        self.assertNotEqual(bin_u1, bin_u50)

        # Tile 0 x-columns [0:8] should have peak at u=1
        tile0_spectrum = fields["s0_u"][0, 0, 0:8, :].sum(axis=0)
        self.assertGreater(tile0_spectrum[bin_u1], 0)
        self.assertEqual(tile0_spectrum[bin_u50], 0)

        # Tile 1 x-columns [8:16] should have peak at u=50
        tile1_spectrum = fields["s0_u"][0, 0, 8:16, :].sum(axis=0)
        self.assertGreater(tile1_spectrum[bin_u50], 0)
        self.assertEqual(tile1_spectrum[bin_u1], 0)

    def test_two_species_distinct_spectra(self):
        """Inject species 0 with |u|=1 and species 1 with |u|=100.
        Verify s0_u and s1_u peak at different bins."""
        P = runko.pic.threeD.ParticleState
        config = make_pic_config(outdir=self.outdir)

        tile_grid = runko.TileGrid(config)
        for idx in tile_grid.local_tile_indices():
            tile = runko.pic.threeD.Tile(idx, config)
            zero = lambda x, y, z: (0, 0, 0)
            tile.set_EBJ(zero, zero, zero)

            def gen0(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(1.0, 0, 0))]
            tile.inject_to_each_cell(0, gen0)

            def gen1(x, y, z):
                return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(100.0, 0, 0))]
            tile.inject_to_each_cell(1, gen1)

            tile_grid.add_tile(tile, idx)
        _ = tile_grid.configure_simulation(config)

        nbins = 100
        umin, umax = 1e-2, 1e4
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=nbins, umin=umin, umax=umax)

        edges = u_bin_edges(hdr)
        bin_u1 = np.searchsorted(edges, 1.0) - 1
        bin_u100 = np.searchsorted(edges, 100.0) - 1

        s0_total = fields["s0_u"].sum(axis=(0, 1, 2))
        s1_total = fields["s1_u"].sum(axis=(0, 1, 2))

        # Species 0 peaks at u=1, not at u=100
        self.assertGreater(s0_total[bin_u1], 0)
        self.assertEqual(s0_total[bin_u100], 0)

        # Species 1 peaks at u=100, not at u=1
        self.assertGreater(s1_total[bin_u100], 0)
        self.assertEqual(s1_total[bin_u1], 0)

    def test_single_species(self):
        """nspecies=1 should produce 4 output fields."""
        config = make_pic_config(outdir=self.outdir, nspecies=1)
        tile_grid = make_grid_with_velocity(config, vx=1.0, nspecies=1)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nspecies=1)

        self.assertEqual(hdr["num_fields"], 4)
        self.assertEqual(hdr["field_names"], ["s0_u", "s0_bx", "s0_by", "s0_bz"])
        self.assertEqual(len(fields), 4)

    def test_three_species(self):
        """nspecies=3 (max for spectra) should produce 12 output fields."""
        config = make_pic_config(outdir=self.outdir, nspecies=3)
        tile_grid = make_grid_with_velocity(config, vx=1.0, nspecies=3)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nspecies=3)

        self.assertEqual(hdr["num_fields"], 12)
        expected_names = []
        for s in range(3):
            expected_names += [f"s{s}_u", f"s{s}_bx", f"s{s}_by", f"s{s}_bz"]
        self.assertEqual(hdr["field_names"], expected_names)
        self.assertEqual(len(fields), 12)

    def test_beta_y_and_z(self):
        """Inject particles with uy-only and uz-only velocities.
        Verify s0_by and s0_bz peaks are correct."""
        config = make_pic_config(outdir=self.outdir)
        nbins = 100

        # Test beta_y: inject uy=1 => gamma=sqrt(2), beta_y=1/sqrt(2)
        tile_grid = make_grid_with_velocity(config, vx=0.0, vy=1.0, vz=0.0)
        hdr, fields = spectra_write_and_read(
            tile_grid, self.outdir, config, nbins=nbins, umin=1e-2, umax=1e2)

        edges = beta_bin_edges(hdr)
        beta_y = 1.0 / np.sqrt(2.0)
        expected_bin = np.searchsorted(edges, beta_y) - 1

        total = fields["s0_by"].sum(axis=(0, 1, 2))
        total_particles = 8**3
        self.assertEqual(total[expected_bin], total_particles)

        # beta_x should peak at 0 (near center bin)
        bx_total = fields["s0_bx"].sum(axis=(0, 1, 2))
        center_bin = nbins // 2
        self.assertEqual(bx_total[center_bin], total_particles)

    def test_header_completeness(self):
        """Verify all header fields."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0)
        hdr, _ = spectra_write_and_read(tile_grid, self.outdir, config)

        self.assertEqual(hdr["magic"], MAGIC)
        self.assertEqual(hdr["version"], 3)
        self.assertEqual(hdr["header_size"], 512)
        self.assertEqual(hdr["dtype_size"], 4)
        self.assertEqual(hdr["stride"], 1)
        self.assertEqual(hdr["lap"], 0)
        self.assertEqual(hdr["Nx"], 1)
        self.assertEqual(hdr["Ny"], 1)
        self.assertEqual(hdr["Nz"], 1)
        self.assertEqual(hdr["NxMesh"], 8)
        self.assertEqual(hdr["NyMesh"], 8)
        self.assertEqual(hdr["NzMesh"], 8)

    def test_nonzero_lap(self):
        """Write with lap=7, verify header contains correct lap."""
        config = make_pic_config(outdir=self.outdir)
        tile_grid = make_grid_with_velocity(config, vx=1.0)
        hdr, _ = spectra_write_and_read(
            tile_grid, self.outdir, config, lap=7)

        self.assertEqual(hdr["lap"], 7)


if __name__ == "__main__":
    unittest.main()
