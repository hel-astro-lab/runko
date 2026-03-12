#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include <mpi.h>

#include "corgi/corgi.h"
#include "core/emf/tile.h"
#include "core/mdgrid_common.h"

namespace mpiio {

/// MPI-IO writer for particle momentum spectra as a function of x.
///
/// Each output file contains a 512-byte header followed by
/// (nspecies * 4) field arrays stored contiguously.  Each field
/// array has shape (Nz, Ny, nx_global, nbins) in C-order where
/// nx_global = Nx * max(1, NxMesh / stride).
///
/// The 4 spectra per species are:
///   u      = sqrt(ux^2 + uy^2 + uz^2)         — log10 bins [umin, umax]
///   beta_x = ux / gamma                        — linear bins [-1, +1]
///   beta_y = uy / gamma                        — linear bins [-1, +1]
///   beta_z = uz / gamma                        — linear bins [-1, +1]
/// where gamma = sqrt(1 + ux^2 + uy^2 + uz^2).
///
/// Per tile, histograms are accumulated using atomic_add on the
/// device tile_buf_, then written via independent MPI_File_write_at
/// (one write per tile per field: nxt*nbins contiguous floats).
template<size_t D>
class SpectraWriter {
  static_assert(D == 3, "Only 3D supported");

public:
  SpectraWriter(
    const std::string& prefix,
    int Nx, int NxMesh,
    int Ny, int NyMesh,
    int Nz, int NzMesh,
    int stride,
    int nbins,
    float umin,
    float umax,
    int nspecies = 2);

  /// Write spectra for all local tiles to prefix/spectra_<lap>.bin
  bool write(corgi::Grid<3>& grid, int lap);

private:
  std::string prefix_;
  int Nx_, Ny_, Nz_;
  int NxMesh_, NyMesh_, NzMesh_;
  int stride_;
  int nbins_;
  float umin_, umax_;
  int nspecies_;
  int num_fields_;      // = nspecies_ * 4
  int nxt_;             // per-tile x output size = max(1, NxMesh / stride)
  int nx_global_;       // = Nx * nxt_

  float log_umin_;      // = log10(umin)
  float inv_dlog_;      // = nbins / (log10(umax) - log10(umin))
  float inv_dbeta_;     // = nbins / 2  (maps [-1,1] to [0, nbins))

  /// Scratch buffer: [nxt_][nbins_] with 12-element SoA per grid point.
  /// Field layout: species * 4 + {0=u, 1=bx, 2=by, 3=bz}.
  runko::SpectraGrid<float> tile_buf_;

  /// Pack (histogram) particle data from one tile into tile_buf_.
  void histogram_tile(emf::Tile<3>& tile);
};

}  // namespace mpiio
