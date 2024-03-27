#pragma once

#include "pic/interpolators/interpolator.h"


namespace pic {

/// Quadratic (2nd order) particle shape interpolator
template<size_t D>
class QuadraticInterpolator :
  public virtual Interpolator<D,3>
{
public: // needs to be public, why is it not public to begin with ?
  void solve(pic::Tile<D>& tile) override;

  double compute( 
        double* /*cx*/, double* /*cy*/, double* /*cz*/, 
        const toolbox::Mesh<float, 3>& /*f*/, 
        const size_t /*iy*/, const size_t /*iz*/,
        int /*i*/, int /*j*/, int /*k*/);

};

} // end of namespace pic




