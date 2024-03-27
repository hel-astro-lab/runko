#pragma once

#include "pic/interpolators/interpolator.h"


namespace pic {

/// 4th order quartic particle shape interpolator
template<size_t D>
class QuarticInterpolator :
  public virtual Interpolator<D,3>
{
  public: 
  void solve(pic::Tile<D>& tile) override;

  double compute( 
        double* /*cx*/, double* /*cy*/, double* /*cz*/, 
        const toolbox::Mesh<float_m, 3>& /*f*/, 
        const size_t /*iy*/, const size_t /*iz*/,
        int /*i*/, int /*j*/, int /*k*/);
  
};


} // end of namespace pic




