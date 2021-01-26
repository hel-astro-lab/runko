#pragma once

#include "propagator.h"

namespace fields {

/// Fourth order staggered finite difference time domain 
// Maxwell's field equation solver.
//
// General formulation with coefficients according to Blinne+2017
//
template<size_t D>
class FDTDGen :
  public virtual Propagator<D>
{
  public:

  /// numerical correction factor to speed of light
  double corr = 1.0; 

  // high-order curl operator coefficients
  //const real_short[6][6][6] = { }
  toolbox::Mesh<real_short,0> CXs;
  toolbox::Mesh<real_short,0> CYs;
  toolbox::Mesh<real_short,0> CZs;

  FDTDGen() :
      CXs(3,3,3),
      CYs(3,3,3),
      CZs(3,3,3)
  {}

  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace fields
