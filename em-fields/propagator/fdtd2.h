#pragma once

#include "propagator.h"

namespace fields {

/// Second order staggered finite difference time domain 
// Maxwell's field equation solver.
template<size_t D>
class FDTD2 :
  public Propagator<D>
{
  public:

  /// numerical correction factor to speed of light
  double corr = 1.0; 

  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace fields
