#pragma once

#include "filter.h"

namespace fields {

/// Digital 2nd order one-pass binomial filter
template<size_t D>
class Binomial2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  //Binomial2(size_t Nx, size_t Ny, size_t Nz) :
  //  Filter<D>(Nx,Ny,Nz) 
  //{}

  void solve(fields::Tile<D>& tile) override;

};



} // end of namespace fields
