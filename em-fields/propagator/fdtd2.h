#pragma once

#include "propagator.h"

namespace fields {

/// Second order staggered finite difference time domain 
// Maxwell's field equation solver.
template<size_t D>
class FDTD2 :
  public virtual Propagator<D>
{
  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace fields
