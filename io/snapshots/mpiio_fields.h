#pragma once

#include <cstddef>
#include <string>

#include "corgi/corgi.h"
#include "core/emf/tile.h"
#include "core/mdgrid_common.h"

namespace mpiio {

/// MPI-IO writer for field snapshots in a flat binary format.
///
/// Each output file contains a 512-byte header followed by num_fields
/// field arrays (ex, ey, ez, bx, by, bz, jx, jy, jz, n0, ..., n{nspecies-1})
/// stored contiguously in row-major order (z outermost, x innermost).
///
/// write() uses independent MPI_File_write_at calls (one per row per field).
/// write_collective() uses MPI_File_set_view + MPI_File_write_all with
/// derived file types, writing directly from the device tile_buf_ without
/// host staging.  On HIP with GPU-aware MPI (OpenMPI 5+), the device
/// pointer is passed straight to MPI_File_write_all.
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
    int stride,
    int nspecies = 2);

  /// Write all local tiles to prefix/flds_<lap>.bin via MPI-IO.
  /// Uses independent row-by-row MPI_File_write_at calls (no temp buffer).
  bool write(corgi::Grid<3>& grid, int lap);

  /// Write all local tiles via collective MPI-IO with derived file types.
  /// Uses MPI_File_set_view + MPI_File_write_all, writing directly from
  /// the device tile_buf_ (no host staging buffer).  On HIP with
  /// GPU-aware MPI the device pointer is passed straight through.
  bool write_collective(corgi::Grid<3>& grid, int lap);

private:
  std::string prefix_;
  int Nx_, Ny_, Nz_;
  int NxMesh_, NyMesh_, NzMesh_;
  int stride_;
  int nspecies_;           // number of particle species (capped at max_species)
  int num_fields_;         // = num_emf_fields + nspecies_
  int nxt_, nyt_, nzt_;    // per-tile output size after stride
  int nx_, ny_, nz_;       // global output size

  /// Iteration driver for for_each_index: (nzt_, nyt_, nxt_) layout_right.
  /// z outermost (dim 0), x innermost (dim 2) matches file format.
  runko::ScalarGrid<float> iter_grid_;

  /// Scratch buffer: [num_fields_ fields][nzt_][nyt_][nxt_] in SoA layout.
  /// Capacity is (num_emf_fields + max_species) = 14 slots per grid point.
  /// On CPU the device buffer is host-accessible; on GPU the device
  /// pointer is passed directly to GPU-aware MPI.
  runko::IOFieldGrid<float> tile_buf_;

  /// Pack tile field data into tile_buf_ using for_each_index kernels.
  /// Data remains on-device; callers are responsible for any D→H sync.
  void pack_tile(emf::Tile<3>& tile);
};

}  // namespace mpiio
