#pragma once

#include "filter.h"

namespace fields {

/// 2D 3-point compensator filter from Birdsall & Langdon
template<size_t D>
class Compensator2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(fields::Tile<D>& tile) override;

};


/// Optimized digital 2nd order one-pass binomial filter
template<size_t D>
class OptBinomial2 :
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

  /// 3-point weight
  double alpha = 0.5;


  void solve(fields::Tile<D>& tile) override;

};


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

  void solve(fields::Tile<D>& tile) override;

};


/// outrolled 2-strided digital 3-point filter
//
template<size_t D>
class Binomial2Strided2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(fields::Tile<D>& tile) override;

};


} // end of namespace fields
