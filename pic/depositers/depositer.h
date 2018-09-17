#pragma once

#include "../tile.h"
#include "../../definitions.h"


namespace pic {

/// General interface for current depositer
template<size_t D, size_t V>
class Depositer
{

  public:

  Depositer() {};

  virtual ~Depositer() = default;

  /// \brief deposit current to grid
  virtual void solve(pic::Tile<D>& ) = 0;

};

} // end of namespace pic











