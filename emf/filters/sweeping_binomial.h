#pragma once

#include "filter.h"

namespace fields {

/// Optimized digital 2nd order one-pass binomial filter
template<size_t D>
class OptBinomial2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(fields::Tile<D>& tile) override;

};


} // end of namespace fields
