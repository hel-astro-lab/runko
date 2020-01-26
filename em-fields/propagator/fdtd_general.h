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
	double coeff1= 9./ 8.0;
	double coeff2=-1./24.0;

  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;
};


} // end of namespace fields
