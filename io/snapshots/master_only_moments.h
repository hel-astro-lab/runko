#pragma once

#include <vector>
#include <string>

#include "io/snapshots/snapshot.h"
#include "external/corgi/corgi.h"


namespace h5io { 

/// IO for calculating particle distribution moments
//
// Calculates full stress tensor, bulk velocites, and number densities
// Uses only the root rank to store the large global array.
template<size_t D>
class MasterPicMomentsWriter :
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

    // global mesh size
    int nx;
    int ny;
    int nz;
      
    // internal mesh tile size
    int nxM;
    int nyM;
    int nzM;

    /// data stride length
    int stride = 1;
      
    // mpi send buffer
    std::vector< toolbox::Mesh<float_m,0> > sbuf;
      
    // status of master rank that is differnet from the rest
    bool master_is_initialized = false;

    /// constructor that creates a name and opens the file handle
    MasterPicMomentsWriter(
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

      // mesh sizes
      nxM = NxMesh/stride;
      nyM = NyMesh/stride;
      nzM = NzMesh/stride;

      nxM = nxM == 0 ? 1 : nxM;
      nyM = nyM == 0 ? 1 : nyM;
      nzM = nzM == 0 ? 1 : nzM;

      // everyone
      sbuf.emplace_back(nxM, nyM, nzM); // only one collective receive buffer
      rbuf.emplace_back(nxM, nyM, nzM); // only one collective receive buffer

    }

    /// read tile meshes into memory
    void read_tiles(corgi::Grid<D>& grid) override;
      
    // local function
    void read_tile_feature(corgi::Grid<D>& grid, uint64_t cid, int ifea);

    /// write hdf5 file
    bool write(corgi::Grid<D>& grid, int lap) override;

    // alternative version of the mpi communication
    void mpi_reduce_snapshots(corgi::Grid<D>& grid) override;
};

} // end of namespace h5io

