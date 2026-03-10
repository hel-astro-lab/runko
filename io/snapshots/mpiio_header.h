#pragma once

#include <bit>
#include <cstdint>
#include <mpi.h>

namespace mpiio {

static_assert(
  std::endian::native == std::endian::little,
  "MPI-IO binary format assumes little-endian");

/// Magic number: "RNKO" in ASCII
inline constexpr uint32_t magic = 0x524E4B4F;

/// Binary format version
inline constexpr uint32_t version = 3;

/// Fixed header size in bytes
inline constexpr int32_t header_size = 512;

/// Maximum number of particle species stored
inline constexpr int32_t max_species = 5;

/// Number of electromagnetic field quantities (ex, ey, ez, bx, by, bz, jx, jy, jz)
inline constexpr int32_t num_emf_fields = 9;

/// EMF field names in canonical order
inline constexpr const char* emf_field_names[9] = {
  "ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz"
};

/// Write the 512-byte binary header to an MPI file handle.
/// Returns MPI error code (MPI_SUCCESS on success).
///
/// Header layout (all little-endian):
///   [0:4]     uint32  magic
///   [4:8]     uint32  version
///   [8:12]    uint32  header_size
///   [12:16]   uint32  num_fields
///   [16:20]   int32   nx  (global output dims after stride)
///   [20:24]   int32   ny
///   [24:28]   int32   nz
///   [28:32]   int32   stride
///   [32:36]   int32   Nx  (tile grid dims)
///   [36:40]   int32   Ny
///   [40:44]   int32   Nz
///   [44:48]   int32   NxMesh
///   [48:52]   int32   NyMesh
///   [52:56]   int32   NzMesh
///   [56:60]   int32   lap
///   [60:64]   uint32  dtype_size = 4
///   [64:64+num_fields*16]  char[16]*num_fields  field names (null-padded)
///   [remainder:512]  reserved zeros
int write_header(
  MPI_File fh,
  int32_t nx, int32_t ny, int32_t nz,
  int32_t stride,
  int32_t Nx, int32_t Ny, int32_t Nz,
  int32_t NxMesh, int32_t NyMesh, int32_t NzMesh,
  int32_t lap,
  int32_t num_fields);

}  // namespace mpiio
