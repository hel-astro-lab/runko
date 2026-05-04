# Copyright 2025 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared constants for runko v5 MPI-IO binary snapshot format."""

MAGIC = 0x524E4B4F          # "RNKO" in little-endian uint32
SUPPORTED_VERSIONS = (3,)
HEADER_SIZE = 512
