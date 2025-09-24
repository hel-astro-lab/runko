#pragma once

#include "definitions.h"
#include "core/emf/tile.h"
#include "tools/mesh.h"


namespace emf {

/// General interface for filters
template<size_t D>
class Filter
{
  public:

  /// grid size along x
  int Nx;
    
  /// grid size along y
  int Ny;

  /// grid size along z
  int Nz;

  ///internal scratch container (size equal to jx/jy/jz)
  toolbox::Mesh<float, 3> tmp;

  Filter(int Nx, int Ny, int Nz) : 
    Nx{Nx}, Ny{Ny}, Nz{Nz},
    tmp{Nx,Ny,Nz}
  {}

  virtual ~Filter() = default;

  virtual void solve(emf::Tile<D>& tile) = 0;

};


} // end of namespace emf


