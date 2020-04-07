#pragma once

#include <vector>
#include <string>

#include "snapshot.h"
#include "../../corgi/corgi.h"


namespace h5io { 

/// IO for calculating particle distribution moments
//
// Calculates full stress tensor, bulk velocites, and number densities
//
template<size_t D>
class PicMomentsWriter :
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
    const string file_name = "moms";

    // internal mesh size
    int nx;
    int ny;
    int nz;

    /// data stride length
    int stride = 1;

    /// constructor that creates a name and opens the file handle
    PicMomentsWriter(
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

      // add correct amount of data containers
      for(size_t i=0; i<14; i++) arrs.emplace_back(nx, ny, nz);
      rbuf.emplace_back(nx, ny, nz);
    }

    /// read tile meshes into memory
    void read_tiles(corgi::Grid<D>& grid) override;

    /// write hdf5 file
    bool write(corgi::Grid<D>& grid, int lap) override;

};

} // end of namespace h5io

