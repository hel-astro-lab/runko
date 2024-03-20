#pragma once

#include "filter.h"

namespace emf {

/// Digital 3-point general filter with variable weight.
//
// Can be made digital filter with alpha = 1/2
// or a compensator for alpha = 3/2
//
template<size_t D>
class General3p :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  /// 3-point weight
  double alpha = 0.5;

  void solve(emf::Tile<D>& tile) override;
};


} // end of namespace emf
