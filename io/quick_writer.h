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
    std::vector< toolbox::Mesh<double> > arrs;


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
      int nx = Nx*NxMesh/stride;
      int ny = Ny*NyMesh/stride;
      int nz = Nz*NzMesh/stride;

      for(size_t i=0; i<10; i++) arrs.emplace_back(nx, ny, nz);
    }

    /// read tile meshes into memory
    void read_tiles(corgi::Node<D>& grid);

    bool write(corgi::Node<D>& grid, int lap);
};


} // end of namespace h5io

