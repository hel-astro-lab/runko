"""Shared constants for runko v5 MPI-IO binary snapshot format."""

MAGIC = 0x524E4B4F          # "RNKO" in little-endian uint32
SUPPORTED_VERSIONS = (3,)
HEADER_SIZE = 512
