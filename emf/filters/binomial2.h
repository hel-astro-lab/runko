#pragma once

#include "filter.h"

namespace emf {

/// Digital 2nd order one-pass binomial filter
template<size_t D>
class Binomial2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(emf::Tile<D>& tile) override;

};

} // end of namespace emf
