#pragma once

#include <cstdint>
#include <cstring>
#include <mpi.h>

namespace mpiio {

/// Magic number: "RNKO" in ASCII
static constexpr uint32_t magic = 0x524E4B4F;

/// Binary format version
static constexpr uint32_t version = 1;

/// Fixed header size in bytes
static constexpr int32_t header_size = 256;

/// Number of field quantities stored
static constexpr int32_t num_fields = 10;

/// Field names in canonical order
static constexpr const char* field_names[10] = {
  "ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz", "rho"
};

/// Write the 256-byte binary header to an MPI file handle.
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
///   [64:224]  char[16]*10  field names (null-padded)
///   [224:256] reserved zeros
inline void write_header(
  MPI_File fh,
  int32_t nx, int32_t ny, int32_t nz,
  int32_t stride,
  int32_t Nx, int32_t Ny, int32_t Nz,
  int32_t NxMesh, int32_t NyMesh, int32_t NzMesh,
  int32_t lap)
{
  char buf[256];
  std::memset(buf, 0, sizeof(buf));

  auto put32 = [&](int offset, uint32_t val) {
    std::memcpy(buf + offset, &val, 4);
  };

  auto puti32 = [&](int offset, int32_t val) {
    std::memcpy(buf + offset, &val, 4);
  };

  put32(0, magic);
  put32(4, version);
  put32(8, static_cast<uint32_t>(header_size));
  put32(12, static_cast<uint32_t>(num_fields));

  puti32(16, nx);
  puti32(20, ny);
  puti32(24, nz);
  puti32(28, stride);
  puti32(32, Nx);
  puti32(36, Ny);
  puti32(40, Nz);
  puti32(44, NxMesh);
  puti32(48, NyMesh);
  puti32(52, NzMesh);
  puti32(56, lap);
  put32(60, 4);  // dtype_size = sizeof(float)

  // Field names: 10 entries of 16 chars each, starting at offset 64
  for (int f = 0; f < num_fields; f++) {
    std::strncpy(buf + 64 + f * 16, field_names[f], 15);
  }

  // [224:256] already zeroed

  MPI_Status status;
  MPI_File_write_at(fh, 0, buf, 256, MPI_BYTE, &status);
}

}  // namespace mpiio
