"""Reader for runko v5 MPI-IO binary spectra snapshot files.

Parses the fixed header and returns numpy arrays of spectral data.

Binary layout
-------------
Bytes 0 to header_size-1 : header (standard RNKO header + spectra metadata)
Bytes header_size+       : num_fields contiguous (nz, ny, nx, nbins) float32 arrays

Spectra metadata in header:
  [256:260]  int32   nbins
  [260:264]  float   umin
  [264:268]  float   umax
Beta limits are always [-1, +1].
"""

import os
import struct
from pathlib import Path

import numpy as np

from runko.mpiio_constants import MAGIC

SUPPORTED_VERSIONS = (3,)


def read_spectra_header(path: str | os.PathLike) -> dict:
    """Parse the header of an MPI-IO spectra snapshot file.

    Returns
    -------
    dict
        Standard header fields plus nbins, umin, umax.
    """
    with open(path, "rb") as fh:
        raw = fh.read(12)
        (magic, ver, hdr_size) = struct.unpack_from("<III", raw, 0)

        if magic != MAGIC:
            raise ValueError(
                f"Bad magic number: expected 0x{MAGIC:08X}, got 0x{magic:08X}"
            )
        if ver not in SUPPORTED_VERSIONS:
            raise ValueError(
                f"Unsupported format version: expected one of {SUPPORTED_VERSIONS}, got {ver}"
            )

        raw += fh.read(hdr_size - 12)

    _hdr_fmt = "<IIIIiiiiiiiiiiII"
    (
        magic, version, header_size, num_fields,
        nx, ny, nz, stride,
        Nx, Ny, Nz,
        NxMesh, NyMesh, NzMesh,
        lap, dtype_size,
    ) = struct.unpack_from(_hdr_fmt, raw, 0)

    field_names = []
    for i in range(num_fields):
        offset = 64 + i * 16
        name_bytes = raw[offset:offset + 16]
        name = name_bytes.split(b"\x00", 1)[0].decode("ascii")
        field_names.append(name)

    # Spectra-specific metadata at byte 256
    (nbins,) = struct.unpack_from("<i", raw, 256)
    (umin,) = struct.unpack_from("<f", raw, 260)
    (umax,) = struct.unpack_from("<f", raw, 264)

    return {
        "magic":       magic,
        "version":     version,
        "header_size": header_size,
        "num_fields":  num_fields,
        "nx":          nx,
        "ny":          ny,
        "nz":          nz,
        "stride":      stride,
        "Nx":          Nx,
        "Ny":          Ny,
        "Nz":          Nz,
        "NxMesh":      NxMesh,
        "NyMesh":      NyMesh,
        "NzMesh":      NzMesh,
        "lap":         lap,
        "dtype_size":  dtype_size,
        "field_names": field_names,
        "nbins":       nbins,
        "umin":        umin,
        "umax":        umax,
    }


def read_spectra_snapshot(path: str | os.PathLike) -> dict[str, np.ndarray]:
    """Read all spectra fields from an MPI-IO spectra snapshot.

    Returns
    -------
    dict[str, numpy.ndarray]
        Mapping from field name to a (nz, ny, nx, nbins) float32 array.
        nz and ny are tile-grid dimensions; sum over axes 0,1 to get
        the total spectrum as a function of x.
    """
    hdr = read_spectra_header(path)
    nx, ny, nz = hdr["nx"], hdr["ny"], hdr["nz"]
    nbins = hdr["nbins"]
    hdr_size = hdr["header_size"]
    field_size = nz * ny * nx * nbins  # float32 elements per field

    fields: dict[str, np.ndarray] = {}
    for idx, name in enumerate(hdr["field_names"]):
        offset = hdr_size + idx * field_size * hdr["dtype_size"]
        arr = np.memmap(
            path,
            dtype=np.float32,
            mode="r",
            offset=offset,
            shape=(nz, ny, nx, nbins),
        )
        fields[name] = np.array(arr)

    return fields


def u_bin_edges(hdr: dict) -> np.ndarray:
    """Return the nbins+1 logarithmically-spaced bin edges for u spectra."""
    return np.logspace(
        np.log10(hdr["umin"]),
        np.log10(hdr["umax"]),
        hdr["nbins"] + 1,
    )


def u_bin_centers(hdr: dict) -> np.ndarray:
    """Return the nbins logarithmic bin centers for u spectra (geometric mean)."""
    edges = u_bin_edges(hdr)
    return np.sqrt(edges[:-1] * edges[1:])


def beta_bin_edges(hdr: dict) -> np.ndarray:
    """Return the nbins+1 linearly-spaced bin edges for beta spectra [-1, +1]."""
    return np.linspace(-1.0, 1.0, hdr["nbins"] + 1)


def beta_bin_centers(hdr: dict) -> np.ndarray:
    """Return the nbins linear bin centers for beta spectra."""
    edges = beta_bin_edges(hdr)
    return 0.5 * (edges[:-1] + edges[1:])
