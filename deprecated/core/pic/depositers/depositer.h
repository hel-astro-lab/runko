#pragma once

#include "core/pic/tile.h"
#include "definitions.h"


namespace pic {

/// General interface for current depositer
template<size_t D, size_t V>
class Depositer
{

  public:

  Depositer() = default;

  virtual ~Depositer() = default;

  /// \brief deposit current to grid
  virtual void solve(pic::Tile<D>& ) = 0;

};

} // end of namespace pic











