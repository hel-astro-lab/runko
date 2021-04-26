#pragma once

#include "interpolator.h"


namespace pic {

/// Linear (1st order) particle shape interpolator
template<size_t D, size_t V>
class LinearInterpolator :
  public virtual Interpolator<D,V>
{
public: // needs to be public, why is it not public to begin with ?
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




