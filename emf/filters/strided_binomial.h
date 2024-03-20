#pragma once

#include "filter.h"

namespace emf {

/// Strided digital 3-point general filter with variable weight.
//
// Can be made digital filter with alpha = 1/2
// or a compensator for alpha = 3/2
//
// Stride parameter allows to damp multiples of Nyquist frequency. 
//
// Ref: Vay, Geddes, Cormier-Michel, Grote, J. Comp. Phys. 2011
//
template<size_t D>
class General3pStrided :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  /// 3-point weight
  double alpha = 0.5;

  /// stride parameter s; s=1 gives standard non-strided
  int stride = 1;

  void solve(emf::Tile<D>& tile) override;

};


/// outrolled 2-strided digital 3-point filter
//
template<size_t D>
class Binomial2Strided2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(emf::Tile<D>& tile) override;

};


} // end of namespace emf
