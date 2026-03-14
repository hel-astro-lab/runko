"""
Multi-rank MPI-IO spectra snapshot writer tests.

Run with: mpirun -np 4 python test_mpiio_spectra.py
"""

import mpi_unittest

import shutil
import tempfile
import numpy as np

from mpi4py import MPI

import runko
from runko_cpp_bindings.emf.threeD import MpiioSpectraWriter
from runko.mpiio_spectra_reader import (
    read_spectra_header, read_spectra_snapshot, u_bin_edges,
)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def broadcast_outdir():
    """Create a temp directory on rank 0 and broadcast its path to all ranks."""
    outdir = None
    if rank == 0:
        outdir = tempfile.mkdtemp(prefix="runko-mpiio-spectra-test-")
    outdir = comm.bcast(outdir, root=0)
    return outdir


def cleanup_outdir(outdir):
    """Remove the temp directory on rank 0."""
    if rank == 0:
        shutil.rmtree(outdir, ignore_errors=True)


def find_spectra_file(outdir):
    """Return path to the spectra output file."""
    import os
    files = [f for f in os.listdir(outdir)
             if f.startswith("pspectra_") and f.endswith(".bin")]
    assert len(files) > 0, "No spectra output file"
    return os.path.join(outdir, files[0])


# ---------------------------------------------------------------------------
# Test 1: constant velocity across all tiles on 4 ranks
# ---------------------------------------------------------------------------
def test_multirank_constant_velocity():
    """2x2x2 grid, 8^3 mesh per tile, hilbert_curve partitioning.
    Every tile has ppc=1 particles with vel=(3,4,0) => |u|=5.
    After writing, rank 0 reads and verifies all particles land in
    the u=5 bin."""

    outdir = broadcast_outdir()
    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nx = 2
    config.Ny = 2
    config.Nz = 2
    config.NxMesh = 8
    config.NyMesh = 8
    config.NzMesh = 8
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.Nt = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.delgam = 1e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"

    P = runko.pic.threeD.ParticleState

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, config)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)

        def gen(x, y, z):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(3.0, 4.0, 0.0))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    nbins = 100
    umin, umax = 1e-2, 1e4
    nspecies = 2
    writer = MpiioSpectraWriter(outdir, 2, 8, 2, 8, 2, 8, 1,
                                nbins, umin, umax, nspecies)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_spectra_file(outdir)
        hdr = read_spectra_header(path)
        fields = read_spectra_snapshot(path)

        # Output shape: (Nz*NzMesh=2, Ny*NyMesh=2, nx=16, nbins)
        assert hdr["nx"] == 16, f"Expected nx=16, got {hdr['nx']}"
        assert hdr["ny"] == 2, f"Expected ny=2, got {hdr['ny']}"
        assert hdr["nz"] == 2, f"Expected nz=2, got {hdr['nz']}"

        edges = u_bin_edges(hdr)
        expected_bin = np.searchsorted(edges, 5.0) - 1

        total_particles = 8 * 8**3  # 8 tiles
        total = fields["s0_u"].sum(axis=(0, 1, 2))

        assert total[expected_bin] == total_particles, \
            f"Expected {total_particles} in bin {expected_bin}, got {total[expected_bin]}"

        # All other bins should be zero
        mask = np.ones(nbins, dtype=bool)
        mask[expected_bin] = False
        np.testing.assert_array_equal(total[mask], 0)

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 2: tile placement — verify each tile's spectra lands at the correct
#          spatial position
# ---------------------------------------------------------------------------
def test_multirank_tile_placement():
    """2x2x1 grid, 8^3 mesh per tile.
    Tile (0,*,*) has |u|=1, tile (1,*,*) has |u|=50.
    Rank 0 reads and verifies x-resolved spectra are correct."""

    outdir = broadcast_outdir()
    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nx = 2
    config.Ny = 2
    config.Nz = 1
    config.NxMesh = 8
    config.NyMesh = 8
    config.NzMesh = 8
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.Nt = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.delgam = 1e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"

    P = runko.pic.threeD.ParticleState

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, config)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)

        vx = 1.0 if idx[0] == 0 else 50.0

        def gen(x, y, z, _vx=vx):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(_vx, 0.0, 0.0))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    nbins = 100
    umin, umax = 1e-2, 1e4
    nspecies = 2
    writer = MpiioSpectraWriter(outdir, 2, 8, 2, 8, 1, 8, 1,
                                nbins, umin, umax, nspecies)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_spectra_file(outdir)
        hdr = read_spectra_header(path)
        fields = read_spectra_snapshot(path)

        edges = u_bin_edges(hdr)
        bin_u1 = np.searchsorted(edges, 1.0) - 1
        bin_u50 = np.searchsorted(edges, 50.0) - 1
        assert bin_u1 != bin_u50

        # Tile x=0 columns [0:8]: peak at u=1
        tile0 = fields["s0_u"][:, :, 0:8, :].sum(axis=(0, 1, 2))
        assert tile0[bin_u1] > 0, f"Tile 0 should have particles at u=1"
        assert tile0[bin_u50] == 0, f"Tile 0 should NOT have particles at u=50"

        # Tile x=1 columns [8:16]: peak at u=50
        tile1 = fields["s0_u"][:, :, 8:16, :].sum(axis=(0, 1, 2))
        assert tile1[bin_u50] > 0, f"Tile 1 should have particles at u=50"
        assert tile1[bin_u1] == 0, f"Tile 1 should NOT have particles at u=1"

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 3: particle count conservation across ranks
# ---------------------------------------------------------------------------
def test_multirank_count_conservation():
    """2x2x2 grid. Verify total particle count in spectra matches
    the expected number across all ranks."""

    outdir = broadcast_outdir()
    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
    config.Nx = 2
    config.Ny = 2
    config.Nz = 2
    config.NxMesh = 8
    config.NyMesh = 8
    config.NzMesh = 8
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.Nt = 1
    config.field_propagator = "FDTD2"
    config.outdir = outdir
    config.q0 = -1
    config.m0 = 1
    config.q1 = 1
    config.m1 = 1
    config.delgam = 1e-5
    config.temperature_ratio = 1.0
    config.sigma = 40
    config.c_omp = 1
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st_atomic"

    P = runko.pic.threeD.ParticleState

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.pic.threeD.Tile(idx, config)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)

        def gen(x, y, z):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(2.0, 3.0, 4.0))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    nbins = 50
    writer = MpiioSpectraWriter(outdir, 2, 8, 2, 8, 2, 8, 1,
                                nbins, 1e-2, 1e4, 2)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_spectra_file(outdir)
        fields = read_spectra_snapshot(path)

        total_particles = 8 * 8**3  # 8 tiles * 512 ppc=1
        u_total = fields["s0_u"].sum()
        assert abs(u_total - total_particles) <= 1, \
            f"Expected {total_particles}, got {u_total}"

        bx_total = fields["s0_bx"].sum()
        assert abs(bx_total - total_particles) <= 1, \
            f"Expected {total_particles}, got {bx_total}"

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


if __name__ == "__main__":
    test_multirank_constant_velocity()
    test_multirank_tile_placement()
    test_multirank_count_conservation()
