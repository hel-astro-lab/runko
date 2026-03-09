"""Shared helpers for MPI-IO field snapshot tests."""

import os

import runko
from runko_cpp_bindings.emf.threeD import MpiioFieldsWriter
from pytools.mpiio_reader import read_header, read_field_snapshot


FIELD_NAMES = ["ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz", "rho"]


def make_config(Nx=1, Ny=1, Nz=1, NxMesh=8, NyMesh=8, NzMesh=8,
                outdir=None, tile_partitioning="catepillar_track"):
    """Create a standard test configuration."""

    config = runko.Configuration(None)
    config.tile_partitioning = tile_partitioning
    if tile_partitioning == "catepillar_track":
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
    """Create a TileGrid, populate tiles with fields, and configure simulation.

    Returns only tile_grid (simulation object is unused by callers).
    """

    tile_grid = runko.TileGrid(config)
    for idx in tile_grid.local_tile_indices():
        tile = runko.emf.threeD.Tile(idx, config)
        tile.set_EBJ(E_func, B_func, J_func)
        tile_grid.add_tile(tile, idx)

    _ = tile_grid.configure_simulation(config)
    return tile_grid


def find_output_file(outdir):
    """Return the path to the first output file in outdir."""

    files = [f for f in os.listdir(outdir)
             if os.path.isfile(os.path.join(outdir, f))]
    assert len(files) > 0, "No output file was created"
    return os.path.join(outdir, files[0])


def write_and_read(tile_grid, outdir, config, stride=1, lap=0):
    """Write a snapshot and read it back via the Python reader.

    Returns (header_dict, fields_dict).
    """

    writer = MpiioFieldsWriter(
        outdir,
        config.Nx, config.NxMesh,
        config.Ny, config.NyMesh,
        config.Nz, config.NzMesh,
        stride)
    writer.write(tile_grid._corgi_grid, lap)
    path = find_output_file(outdir)
    hdr = read_header(path)
    fields = read_field_snapshot(path)
    return hdr, fields
