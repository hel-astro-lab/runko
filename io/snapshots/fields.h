#pragma once

#include <vector>
#include <string>

#include "snapshot.h"
#include "../../corgi/corgi.h"


namespace h5io { 


template<size_t D>
class FieldsWriter :
  public SnapshotWriter<D>
{

  public:

    using SnapshotWriter<D>::fname;
    using SnapshotWriter<D>::extension;
    using SnapshotWriter<D>::arrs;
    using SnapshotWriter<D>::rbuf;

    using SnapshotWriter<D>::mpi_reduce_snapshots;

  public:

    /// general file name used for outputs
    const string file_name = "flds";

    // internal mesh size
    int nx;
    int ny;
    int nz;

    /// data stride length
    int stride = 1;

    /// constructor that creates a name and opens the file handle
    FieldsWriter(
        const std::string& prefix, 
        int Nx, int NxMesh,
        int Ny, int NyMesh,
        int Nz, int NzMesh,
        int stride) :
      SnapshotWriter<D>{prefix},
      stride{stride}
    {

      //fname = prefix + "-" + to_string(lap) + extension;
      nx = Nx*NxMesh/stride;
      ny = Ny*NyMesh/stride;
      nz = Nz*NzMesh/stride;

      nx = nx == 0 ? 1 : nx;
      ny = ny == 0 ? 1 : ny;
      nz = nz == 0 ? 1 : nz;

      for(size_t i=0; i<10; i++) arrs.emplace_back(nx, ny, nz);
      for(size_t i=0; i<10; i++) rbuf.emplace_back(nx, ny, nz);
    }

    /// read tile meshes into memory
    void read_tiles(corgi::Grid<D>& grid) override;

    /// write hdf5 file
    bool write(corgi::Grid<D>& grid, int lap) override;

};

} // end of namespace h5io

