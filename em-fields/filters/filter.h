#pragma once

#include "../tile.h"
#include "../../definitions.h"


namespace fields {

/// General interface for filters
template<size_t D>
class Filter
{
  public:

  Filter() {};

  virtual ~Filter() = default;

  virtual void solve(Tile<D>& tile) = 0;

};


} // end of namespace fields


