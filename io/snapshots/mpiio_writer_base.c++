#include "io/snapshots/mpiio_writer_base.h"

#include <mpi.h>


template<>
bool mpiio::WriterBase<3>::write(corgi::Grid<3>& grid, int lap)
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


//--------------------------------------------------
// explicit template class instantiation
template class mpiio::WriterBase<3>;
