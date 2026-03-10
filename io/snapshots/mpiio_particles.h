#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include <mpi.h>

#include "corgi/corgi.h"
#include "core/mdgrid_common.h"
#include "io/snapshots/mpiio_header.h"

namespace mpiio {

/// MPI-IO writer for sampled particle snapshots in a flat binary format.
///
/// Each output file contains a 512-byte header followed by num_prtcl_fields
/// 1D arrays (x, y, z, ux, uy, uz, ex, ey, ez, bx, by, bz) stored
/// contiguously in SoA layout.
///
/// Writes a fair sub-sample of n_prtcls particles (user-specified total)
/// with per-tile quotas proportional to each tile's particle count.
/// Particles do not have unique IDs; the sample changes every snapshot.
///
/// Uses independent MPI_File_write_at calls (one per field per rank).
/// Each rank computes its file offset via MPI_Exscan on the local
/// sampled count.
///
/// Sampling, position-to-global conversion, and field interpolation
/// are computed on-device using tyvi for_each_index kernels.
template<size_t D>
class ParticlesWriter {
  static_assert(D == 3, "Only 3D supported");

public:
  ParticlesWriter(const std::string& prefix, int64_t n_prtcls, int species = 0);

  /// Write sampled particles to prefix/prtcls_<species>_<lap>.bin via MPI-IO.
  bool write(corgi::Grid<3>& grid, int lap);

private:
  std::string prefix_;
  int64_t n_prtcls_;               // target total sample count
  int species_;                    // which particle species to sample

  static constexpr int num_prtcl_fields_ = mpiio::num_prtcl_fields;

  /// Scratch buffer: [local_sampled][12] in SoA layout.
  /// Field f is contiguous at offset f * local_sampled in the underlying span.
  runko::PrtclFieldList<float> prtcl_buf_;
};

}  // namespace mpiio
