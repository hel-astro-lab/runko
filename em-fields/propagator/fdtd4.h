#pragma once

#include "propagator.h"

namespace fields {

/// Fourth order staggered finite difference time domain 
// Maxwell's field equation solver.
template<size_t D>
class FDTD4 :
  public Propagator<D>
{
  public:

  /// numerical correction factor to speed of light
  double corr = 1.0; 

  // high-order curl operator coefficients
	double coeff1= 9./ 8.0;
	double coeff2=-1./24.0;

  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace fields
