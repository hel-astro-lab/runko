#pragma once

#include "propagator.h"
#include "../../ffe/tile.h"

namespace fields {

/// Second order staggered finite difference time domain 
// Maxwell's field equation solver.
template<size_t D>
class FDTD2_pml :
  public Propagator<D>
{
  public:

  /// numerical correction factor to speed of light
  double corr = 1.0; 

  /// global grid center
  double cenx, ceny, cenz;

  /// grid dimension scales
  double radx, rady, radz;

  /// reference damping radius (in units of normalized unit sphere)
  double rad_lim = 0.9;

  /// absorption coefficient (should be cfl/3)
  double norm_abs = 0.1;

  /// damping function
  real_short lambda(real_short sx, real_short sy, real_short sz);



  void push_e(Tile<D>& tile) override;

  void push_half_b(Tile<D>& tile) override;

  /// FFE simultaneous E & B update with PML damping
  void push_eb(::ffe::Tile<D>& tile);
};


} // end of namespace fields
