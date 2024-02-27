#pragma once

#include <vector>
#include <string>

#include "snapshot.h"
#include "../../corgi/corgi.h"


namespace h5io { 


/// IO object for storing peripheral of a 3D simulation domain
class FieldSliceWriter :
  public SnapshotWriter<3>
{
  public:

    using SnapshotWriter<3>::fname;
    using SnapshotWriter<3>::extension;
    using SnapshotWriter<3>::arrs;
    using SnapshotWriter<3>::rbuf;
    using SnapshotWriter<3>::mpi_reduce_snapshots;

  public:

    /// general file name used for outputs
    const string file_name = "slices";

    // internal mesh size
    int nx;
    int ny;
    int nz;

    /// data stride length
    int stride = 1;

    /// slice orientation
    int mode = 0;
       
    /// slice index
    int ind = 0;

    /// constructor that creates a name and opens the file handle
    FieldSliceWriter(
        const std::string& prefix, 
        int Nx, int NxMesh,
        int Ny, int NyMesh,
        int Nz, int NzMesh,
        int n_slices, int mode,
        int stride) :
      SnapshotWriter<3>{prefix},
      stride{stride},
      mode{mode}
    {

      // x-y plane, z = 0
      if(mode == 0) {
        nx = Nx*NxMesh/stride;
        ny = Ny*NyMesh/stride;

      // x-z plane, y = 0
      } else if(mode == 1) {

        nx = Nx*NxMesh/stride;
        ny = Nz*NzMesh/stride;

      // y-z plane, x = 0
      } else if(mode == 2) {
        nx = Ny*NyMesh/stride;
        ny = Nz*NzMesh/stride;
      }

      nz = n_slices;

      nx = nx == 0 ? 1 : nx;
      ny = ny == 0 ? 1 : ny;

      for(size_t i=0; i<10; i++) arrs.emplace_back(nx, ny, nz);
      rbuf.emplace_back(nx, ny, nz); // only one collective receive buffer
    }

    /// read tile meshes into memory
    void read_tiles(corgi::Grid<3>& grid) override;

    /// write hdf5 file
    bool write(corgi::Grid<3>& grid, int lap) override;

};

} // end of namespace h5io

