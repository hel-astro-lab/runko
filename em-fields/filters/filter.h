#pragma once

#include "../tile.h"
#include "../../tools/mesh.h"
#include "../../definitions.h"


namespace fields {

/// General interface for filters
template<size_t D>
class Filter
{
  public:

  /// grid size along x
  size_t Nx;
    
  /// grid size along y
  size_t Ny;

  /// grid size along z
  size_t Nz;

  ///internal scratch container (size equal to jx/jy/jz)
  toolbox::Mesh<float_t, 3> tmp;

  Filter(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    tmp(Nx,Ny,Nz)
  {}

  virtual ~Filter() = default;

  virtual void solve(fields::Tile<D>& tile) = 0;

};


} // end of namespace fields


