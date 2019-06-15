#pragma once

#include "../tile.h"
#include "../../definitions.h"


namespace ffe {

/// General interface for FFE current calculations
template<size_t D>
class Current
{
  public:

  Current() {};

  virtual ~Current() = default;

  virtual void solve(Tile<D>& tile) = 0;

};


} // end of namespace ffe

