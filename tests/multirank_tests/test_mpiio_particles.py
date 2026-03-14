"""
Multi-rank MPI-IO particle snapshot writer tests.

Run with: mpirun -np 4 python test_mpiio_particles.py
"""

import mpi_unittest

import shutil
import tempfile
import numpy as np

from mpi4py import MPI

import runko
from runko_cpp_bindings.emf.threeD import MpiioParticlesWriter
from runko.mpiio_prtcl_reader import read_prtcl_header, read_prtcl_snapshot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def broadcast_outdir():
    """Create a temp directory on rank 0 and broadcast its path to all ranks."""
    outdir = None
    if rank == 0:
        outdir = tempfile.mkdtemp(prefix="runko-mpiio-prtcl-test-")
    outdir = comm.bcast(outdir, root=0)
    return outdir


def cleanup_outdir(outdir):
    """Remove the temp directory on rank 0."""
    if rank == 0:
        shutil.rmtree(outdir, ignore_errors=True)


def find_prtcl_file(outdir, species=0):
    """Return path to the particle snapshot file for given species."""
    import os
    files = [f for f in os.listdir(outdir)
             if f.startswith(f"prtcls_{species}_") and f.endswith(".bin")]
    assert len(files) > 0, f"No particle output file for species {species}"
    return os.path.join(outdir, files[0])


# ---------------------------------------------------------------------------
# Test 1: constant velocity across all tiles on 4 ranks
# ---------------------------------------------------------------------------
def test_multirank_constant_velocity():
    """2x2x2 grid, 8^3 mesh per tile, hilbert_curve partitioning.
    Every tile has ppc=1 particles with vel=(1,2,3).
    After writing, rank 0 reads and verifies all particle velocities."""

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
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(1.0, 2.0, 3.0))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    total = 8 * 8**3  # 8 tiles * 512 particles each
    writer = MpiioParticlesWriter(outdir, total * 10, 0)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_prtcl_file(outdir)
        hdr = read_prtcl_header(path)
        fields = read_prtcl_snapshot(path)

        assert hdr["n_prtcls"] == total, \
            f"Expected {total} particles, got {hdr['n_prtcls']}"

        np.testing.assert_allclose(fields["ux"], 1.0, atol=1e-5)
        np.testing.assert_allclose(fields["uy"], 2.0, atol=1e-5)
        np.testing.assert_allclose(fields["uz"], 3.0, atol=1e-5)

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 2: tile placement — verify particles from each tile have correct
#          global positions when written from different ranks
# ---------------------------------------------------------------------------
def test_multirank_tile_placement():
    """2x2x2 grid, 8^3 mesh per tile.
    Each tile injects particles with ux = i+1, uy = j+1, uz = k+1.
    Rank 0 reads and verifies particle positions and velocities."""

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

        vx = float(idx[0] + 1)
        vy = float(idx[1] + 1)
        vz = float(idx[2] + 1)

        def gen(x, y, z, _vx=vx, _vy=vy, _vz=vz):
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(_vx, _vy, _vz))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    total = 8 * 8**3
    writer = MpiioParticlesWriter(outdir, total * 10, 0)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_prtcl_file(outdir)
        hdr = read_prtcl_header(path)
        fields = read_prtcl_snapshot(path)

        assert hdr["n_prtcls"] == total

        # All positions should be non-negative
        assert np.all(fields["x"] >= 0), "Negative x positions"
        assert np.all(fields["y"] >= 0), "Negative y positions"
        assert np.all(fields["z"] >= 0), "Negative z positions"

        # Use velocity encoding to identify tiles and verify each has
        # the expected particle count and distinct spatial positions.
        # Velocity code: ux = i+1, uy = j+1, uz = k+1
        velocity_codes = set()
        for i_tile in range(2):
            for j_tile in range(2):
                for k_tile in range(2):
                    expected_ux = float(i_tile + 1)
                    expected_uy = float(j_tile + 1)
                    expected_uz = float(k_tile + 1)
                    in_tile = (
                        (np.abs(fields["ux"] - expected_ux) < 1e-4) &
                        (np.abs(fields["uy"] - expected_uy) < 1e-4) &
                        (np.abs(fields["uz"] - expected_uz) < 1e-4)
                    )
                    n_in = np.sum(in_tile)
                    assert n_in == 8**3, \
                        f"Tile ({i_tile},{j_tile},{k_tile}): expected {8**3}, got {n_in}"

                    # Each tile's particles should span exactly NxMesh cells
                    x_range = fields["x"][in_tile].max() - fields["x"][in_tile].min()
                    assert abs(x_range - 7.0) < 0.1, \
                        f"Tile ({i_tile},{j_tile},{k_tile}) x range={x_range}"

                    velocity_codes.add((expected_ux, expected_uy, expected_uz))

        assert len(velocity_codes) == 8, \
            f"Expected 8 distinct tiles, got {len(velocity_codes)}"

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


# ---------------------------------------------------------------------------
# Test 3: asymmetric mesh dimensions across multiple ranks
# ---------------------------------------------------------------------------
def test_multirank_asymmetric_mesh():
    """2x2x2 grid, 10x11x13 mesh per tile.
    Constant velocity. Rank 0 verifies total particle count and positions."""

    outdir = broadcast_outdir()
    config = runko.Configuration(None)
    config.tile_partitioning = "hilbert_curve"
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
            return [P(pos=(x + 0.5, y + 0.5, z + 0.5), vel=(5.0, 6.0, 7.0))]
        for species in range(2):
            tile.inject_to_each_cell(species, gen)
        tile_grid.add_tile(tile, idx)
    _ = tile_grid.configure_simulation(config)

    total = 8 * 10 * 11 * 13  # 8 tiles
    writer = MpiioParticlesWriter(outdir, total * 10, 0)
    writer.write(tile_grid._corgi_grid, 0)

    if rank == 0:
        path = find_prtcl_file(outdir)
        hdr = read_prtcl_header(path)
        fields = read_prtcl_snapshot(path)

        assert hdr["n_prtcls"] == total, \
            f"Expected {total}, got {hdr['n_prtcls']}"

        # Positions should be non-negative
        assert np.all(fields["x"] >= 0), "Negative x positions"
        assert np.all(fields["y"] >= 0), "Negative y positions"
        assert np.all(fields["z"] >= 0), "Negative z positions"

        np.testing.assert_allclose(fields["ux"], 5.0, atol=1e-5)
        np.testing.assert_allclose(fields["uy"], 6.0, atol=1e-5)
        np.testing.assert_allclose(fields["uz"], 7.0, atol=1e-5)

    comm.Barrier()
    mpi_unittest.assertEqual(True, True)

    cleanup_outdir(outdir)


if __name__ == "__main__":
    test_multirank_constant_velocity()
    test_multirank_tile_placement()
    test_multirank_asymmetric_mesh()
