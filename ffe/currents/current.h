#pragma once

#include "../tile.h"
#include "../../definitions.h"


namespace ffe {

/// General interface for FFE current calculations
template<size_t D>
class Current
{
  public:

  int Nx;
  int Ny;
  int Nz;

  /// half-cell staggered B
  toolbox::Mesh<double, 1> bxf;
  toolbox::Mesh<double, 1> byf;
  toolbox::Mesh<double, 1> bzf;


  Current(int Nx, int Ny, int Nz) :
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz)
  {};

  virtual ~Current() = default;

  virtual void comp_drift_cur(Tile<D>& tile) = 0;

  virtual void comp_parallel_cur(Tile<D>& tile) = 0;

  virtual void limiter(Tile<D>& tile) = 0;

};


} // end of namespace ffe

