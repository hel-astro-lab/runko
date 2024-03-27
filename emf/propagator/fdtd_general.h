#pragma once

#include "emf/propagator/propagator.h"

namespace emf {

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
  //const float_m[6][6][6] = { }
  toolbox::Mesh<float_m,0> CXs;
  toolbox::Mesh<float_m,0> CYs;
  toolbox::Mesh<float_m,0> CZs;

  FDTDGen() :
      CXs(4,4,4),
      CYs(4,4,4),
      CZs(4,4,4)
  {}

  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace emf
