#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"


namespace pic {

/// General interface for particle interpolators
/// Interpolates electromagnetic fields to particle locations
template<size_t D, size_t V>
class Interpolator
{

  public:

  Interpolator() = default;

  virtual ~Interpolator() = default;

  /// \brief interpolate electromagnetic fields to particle locations
  virtual void solve(pic::Tile<D>& ) = 0;

};


} // end of namespace pic
