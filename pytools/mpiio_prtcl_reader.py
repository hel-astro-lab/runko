"""Reader for runko v5 MPI-IO binary particle snapshot files.

Parses the fixed header and provides zero-copy access to
particle data via numpy memory mapping.

Binary layout
-------------
Bytes 0 to header_size-1 : header (see read_prtcl_header for definitions)
Bytes header_size+       : num_fields contiguous 1D float32 arrays of length n_prtcls
"""

import os
import struct

import numpy as np


# -- module-level constants ---------------------------------------------------

from pytools.mpiio_constants import MAGIC, SUPPORTED_VERSIONS


# -- public API ---------------------------------------------------------------

def read_prtcl_header(path: str | os.PathLike) -> dict:
    """Parse the header of an MPI-IO particle snapshot file.

    Parameters
    ----------
    path : str | os.PathLike
        Path to the binary snapshot file.

    Returns
    -------
    dict
        Header fields: magic, version, header_size, num_fields,
        n_prtcls, species, lap, dtype_size, field_names.

    Raises
    ------
    ValueError
        If the magic number or format version does not match.
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

    (num_fields,) = struct.unpack_from("<I", raw, 12)

    # n_prtcls stored as int64 at offset 16 (occupies the nx+ny slots)
    (n_prtcls,) = struct.unpack_from("<q", raw, 16)

    # species stored in the stride slot at offset 28
    (species,) = struct.unpack_from("<i", raw, 28)

    (lap,) = struct.unpack_from("<i", raw, 56)
    (dtype_size,) = struct.unpack_from("<I", raw, 60)

    field_names = []
    for i in range(num_fields):
        offset = 64 + i * 16
        name_bytes = raw[offset:offset + 16]
        name = name_bytes.split(b"\x00", 1)[0].decode("ascii")
        field_names.append(name)

    return {
        "magic":       magic,
        "version":     ver,
        "header_size": hdr_size,
        "num_fields":  num_fields,
        "n_prtcls":    n_prtcls,
        "species":     species,
        "lap":         lap,
        "dtype_size":  dtype_size,
        "field_names": field_names,
    }


def read_prtcl_snapshot(path: str | os.PathLike) -> dict[str, np.ndarray]:
    """Read all particle fields from an MPI-IO snapshot via memory mapping.

    Parameters
    ----------
    path : str | os.PathLike
        Path to the binary particle snapshot file.

    Returns
    -------
    dict[str, numpy.ndarray]
        Mapping from field name to a 1D float32 array of length n_prtcls.
        Arrays are backed by a read-only memory map (zero-copy).
    """
    hdr = read_prtcl_header(path)
    n = hdr["n_prtcls"]
    hdr_size = hdr["header_size"]

    fields: dict[str, np.ndarray] = {}
    for idx, name in enumerate(hdr["field_names"]):
        offset = hdr_size + idx * n * hdr["dtype_size"]
        arr = np.memmap(
            path,
            dtype=np.float32,
            mode="r",
            offset=offset,
            shape=(n,),
        )
        fields[name] = arr

    return fields
