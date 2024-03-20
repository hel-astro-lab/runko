#pragma once

#include "filter.h"

namespace emf {

/// 2D 3-point compensator filter from Birdsall & Langdon
template<size_t D>
class Compensator2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(emf::Tile<D>& tile) override;

};

} // end of namespace emf
