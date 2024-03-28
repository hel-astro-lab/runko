#pragma once

#include "core/pic/tile.h"
#include "definitions.h"
#include "external/iter/allocator.h"

namespace pic {

/// General interface for particle interpolators
/// Interpolates electromagnetic fields to particle locations
template<size_t D, size_t V>
class Interpolator: public ManagedParent
{

  public:

  Interpolator() = default;

  virtual ~Interpolator() = default;

  /// \brief interpolate electromagnetic fields to particle locations
  virtual void solve(pic::Tile<D>& ) = 0;

};


} // end of namespace pic
