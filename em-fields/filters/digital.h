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

  void solve(fields::Tile<D>& tile) override;

};


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

  double alpha = 0.5;

  void solve(fields::Tile<D>& tile) override;

};



} // end of namespace fields
