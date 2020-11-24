#pragma once

#include "../tile.h"
#include "../../tools/mesh.h"
#include "../../definitions.h"
#include "../../tools/iter/managed_alloc.h"


namespace fields {

/// General interface for filters
template<size_t D>
class Filter: public ManagedParent
{
  public:

  /// grid size along x
  int Nx;
    
  /// grid size along y
  int Ny;

  /// grid size along z
  int Nz;

  ///internal scratch container (size equal to jx/jy/jz)
  toolbox::Mesh<real_short, 3> tmp;

  Filter(int Nx, int Ny, int Nz) : 
    Nx{Nx}, Ny{Ny}, Nz{Nz},
    tmp{Nx,Ny,Nz}
  {}

  virtual ~Filter() = default;

  virtual void solve(fields::Tile<D>& tile) = 0;

};


} // end of namespace fields


