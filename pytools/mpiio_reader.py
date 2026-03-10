"""Reader for runko v5 MPI-IO binary field snapshot files.

Parses the fixed header and provides zero-copy access to
field data via numpy memory mapping.

Binary layout
-------------
Bytes 0 to header_size-1 : header (see read_header for field definitions)
Bytes header_size+       : num_fields contiguous (nz, ny, nx) float32 arrays in C-order
"""

import os
import struct
from pathlib import Path

import numpy as np


# -- module-level constants ---------------------------------------------------

from pytools.mpiio_constants import MAGIC

SUPPORTED_VERSIONS = (1, 2, 3)  # field reader supports older format versions too


# -- public API ---------------------------------------------------------------

def read_header(path: str | os.PathLike) -> dict:
    """Parse the header of an MPI-IO snapshot file.

    Parameters
    ----------
    path : str | os.PathLike
        Path to the binary snapshot file.

    Returns
    -------
    dict
        Header fields: magic, version, header_size, num_fields,
        nx, ny, nz, stride, Nx, Ny, Nz, NxMesh, NyMesh, NzMesh,
        lap, dtype_size, field_names.

    Raises
    ------
    ValueError
        If the magic number or format version does not match.
    """
    # Read just the header bytes (not the whole file, which can be huge).
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

        # Read the remainder of the header
        raw += fh.read(hdr_size - 12)

    # 16 x 32-bit values = 64 bytes:
    #   I  magic       I  version     I  header_size  I  num_fields
    #   i  nx          i  ny          i  nz           i  stride
    #   i  Nx          i  Ny          i  Nz
    #   i  NxMesh      i  NyMesh      i  NzMesh
    #   i  lap         I  dtype_size
    _hdr_fmt = "<IIIIiiiiiiiiiiII"
    assert struct.calcsize(_hdr_fmt) == 64

    (
        magic, version, header_size, num_fields,
        nx, ny, nz, stride,
        Nx, Ny, Nz,
        NxMesh, NyMesh, NzMesh,
        lap, dtype_size,
    ) = struct.unpack_from(_hdr_fmt, raw, 0)

    # field names: num_fields x 16-byte null-padded ASCII strings at offset 64
    field_names = []
    for i in range(num_fields):
        offset = 64 + i * 16
        name_bytes = raw[offset:offset + 16]
        name = name_bytes.split(b"\x00", 1)[0].decode("ascii")
        field_names.append(name)

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
    }


def read_field_snapshot(path: str | os.PathLike) -> dict[str, np.ndarray]:
    """Read all fields from an MPI-IO snapshot via memory mapping.

    Parameters
    ----------
    path : str | os.PathLike
        Path to the binary snapshot file.

    Returns
    -------
    dict[str, numpy.ndarray]
        Mapping from field name to a (nz, ny, nx) float32 array.
        Arrays are backed by a read-only memory map (zero-copy).
    """
    hdr = read_header(path)
    nx, ny, nz = hdr["nx"], hdr["ny"], hdr["nz"]
    hdr_size = hdr["header_size"]
    field_size = nx * ny * nz  # number of float32 elements per field

    fields: dict[str, np.ndarray] = {}
    for idx, name in enumerate(hdr["field_names"]):
        offset = hdr_size + idx * field_size * hdr["dtype_size"]
        arr = np.memmap(
            path,
            dtype=np.float32,
            mode="r",
            offset=offset,
            shape=(nz, ny, nx),
        )
        fields[name] = arr

    return fields


def read_field(path: str | os.PathLike, field_name: str) -> np.ndarray:
    """Read a single named field without loading all fields.

    Parameters
    ----------
    path : str | os.PathLike
        Path to the binary snapshot file.
    field_name : str
        Name of the field to read (must match a header field name).

    Returns
    -------
    numpy.ndarray
        A (nz, ny, nx) float32 array (read-only memory map).

    Raises
    ------
    KeyError
        If *field_name* is not found among the header field names.
    """
    hdr = read_header(path)

    try:
        idx = hdr["field_names"].index(field_name)
    except ValueError:
        raise KeyError(
            f"Field '{field_name}' not found. "
            f"Available fields: {hdr['field_names']}"
        ) from None

    nx, ny, nz = hdr["nx"], hdr["ny"], hdr["nz"]
    hdr_size = hdr["header_size"]
    field_size = nx * ny * nz
    offset = hdr_size + idx * field_size * hdr["dtype_size"]

    return np.memmap(
        path,
        dtype=np.float32,
        mode="r",
        offset=offset,
        shape=(nz, ny, nx),
    )
