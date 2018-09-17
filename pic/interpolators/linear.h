#pragma once

#include "interpolator.h"


namespace pic {

/// Linear (1st order) particle shape interpolator
template<size_t D, size_t V>
class LinearInterpolator :
  public virtual Interpolator<D,V>
{

  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




