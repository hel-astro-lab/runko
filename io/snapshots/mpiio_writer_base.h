#pragma once

#include <cstddef>
#include <string>

#include <mpi.h>

#include "corgi/corgi.h"

namespace mpiio {

/// Base class for MPI-IO snapshot writers (NVI pattern).
///
/// The public write() method implements the shared skeleton:
///   1. prepare_()          — optional pre-write MPI collectives
///   2. make_filename_()    — build output path
///   3. MPI_File_open       — collective open
///   4. write_header_()     — rank 0 writes 512-byte header
///   5. write_payload_()    — all ranks pre-allocate and write their data
///   6. MPI_File_close      — collective close
///
/// Derived classes override the private virtual hooks.
template<size_t D>
class WriterBase {
  static_assert(D == 3, "Only 3D supported");

public:
  WriterBase(const std::string& prefix, int num_fields)
    : prefix_(prefix), num_fields_(num_fields) {}

  virtual ~WriterBase() = default;

  /// Main entry point: prepare -> open -> header -> data -> close.
  bool write(corgi::Grid<3>& grid, int lap)
  {
    // 1. Optional pre-write communication
    if (!prepare_(grid, lap)) return false;

    // 2. Build filename
    std::string filename = make_filename_(lap);

    // 3. Open file (collective)
    MPI_File fh;
    int rc = MPI_File_open(
      MPI_Comm(grid.comm), filename.c_str(),
      MPI_MODE_CREATE | MPI_MODE_WRONLY,
      MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) return false;

    // 4. Rank 0 writes header
    if (grid.comm.rank() == 0) {
      rc = write_header_(fh, lap);
      if (rc != MPI_SUCCESS) { MPI_File_close(&fh); return false; }
    }

    // 5. All ranks write their data
    bool ok = write_payload_(fh, grid);

    // 6. Close
    rc = MPI_File_close(&fh);
    return ok && (rc == MPI_SUCCESS);
  }

protected:
  std::string prefix_;
  int num_fields_;

private:
  /// Build the output filename for this lap.
  virtual std::string make_filename_(int lap) const = 0;

  /// Optional pre-write setup (e.g., MPI collectives for particle counting).
  /// Called before MPI_File_open. Return false to abort.
  virtual bool prepare_(corgi::Grid<3>& /*grid*/, int /*lap*/) { return true; }

  /// Write the 512-byte header. Called on rank 0 only.
  virtual int write_header_(MPI_File fh, int lap) = 0;

  /// Write all per-rank payload data. File is open, header written.
  virtual bool write_payload_(MPI_File fh, corgi::Grid<3>& grid) = 0;
};

}  // namespace mpiio
