#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "corgi/corgi.h"
#include "core/emf/tile.h"

namespace mpiio {

/// MPI-IO writer for field snapshots in a flat binary format.
///
/// Each output file contains a 256-byte header followed by 10 field
/// arrays (ex, ey, ez, bx, by, bz, jx, jy, jz, rho) stored
/// contiguously in row-major order (z outermost, x innermost).
///
/// Every MPI rank writes its own tiles directly at the correct file
/// offset using independent MPI_File_write_at calls, avoiding the
/// reduce-to-rank-0 bottleneck of the HDF5 writer.
template<size_t D>
class FieldsWriter {
  static_assert(D == 3, "Only 3D supported");

public:
  FieldsWriter(
    const std::string& prefix,
    int Nx, int NxMesh,
    int Ny, int NyMesh,
    int Nz, int NzMesh,
    int stride);

  /// Write all local tiles to prefix/flds_<lap>.bin via MPI-IO.
  bool write(corgi::Grid<3>& grid, int lap);

private:
  std::string prefix_;
  int Nx_, Ny_, Nz_;
  int NxMesh_, NyMesh_, NzMesh_;
  int stride_;
  int nxt_, nyt_, nzt_;  // per-tile output size after stride
  int nx_, ny_, nz_;     // global output size

  /// Scratch buffer: [num_fields][nzt_][nyt_][nxt_] row-major
  std::vector<float> tile_buf_;

  /// Fill tile_buf_ from a single emf tile's Yee lattice data.
  void read_tile(emf::Tile<3>& tile);
};

}  // namespace mpiio
