#pragma once

#include <cstddef>
#include <string>

#include "corgi/corgi.h"
#include "core/emf/tile.h"
#include "core/mdgrid_common.h"

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
///
/// Stride subsampling and J volume-sum are computed on-device using
/// tyvi for_each_index kernels (SIMD on CPU, GPU kernels on HIP).
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

  /// Iteration driver for for_each_index: (nzt_, nyt_, nxt_) layout_right.
  /// z outermost (dim 0), x innermost (dim 2) matches file format.
  runko::ScalarGrid<float> iter_grid_;

  /// Scratch buffer: [10 fields][nzt_][nyt_][nxt_] in SoA layout.
  /// On CPU the device buffer is host-accessible; on GPU sync_to_staging
  /// copies to the internal staging buffer for MPI writes.
  runko::IOFieldGrid<float> tile_buf_;

  /// Pack tile field data into tile_buf_ using for_each_index kernels.
  void pack_tile(emf::Tile<3>& tile);
};

}  // namespace mpiio
