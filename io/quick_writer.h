#pragma once

#include <vector>
#include <string>

#include "namer.h"
#include "../corgi/corgi.h"
#include "../tools/mesh.h"


namespace h5io { 


template<size_t D>
class QuickWriter {

  private:

    /// general file extension to be appended to file names
    const string extension = ".h5";
    
    /// Object to handle file names and extensions
    std::string fname;

    /// meshes
    std::vector< toolbox::Mesh<float> > arrs;

    /// mpi receive buffer
    std::vector< toolbox::Mesh<float> > rbuf;

    /// internal image size
    int nx,ny,nz;


  public:

    /// data stride length
    int stride = 1;

    /// constructor that creates a name and opens the file handle
    QuickWriter(
        const std::string& prefix, 
        int Nx, int NxMesh,
        int Ny, int NyMesh,
        int Nz, int NzMesh,
        int stride) :
      fname(prefix),
      stride(stride)
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
    void read_tiles(corgi::Grid<D>& grid);

    /// communicate snapshots with a B-tree cascade to rank 0
    void mpi_reduce_snapshots(corgi::Grid<D>& grid);

    /// write hdf5 file
    bool write(corgi::Grid<D>& grid, int lap);
};


} // end of namespace h5io

