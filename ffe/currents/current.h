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

  virtual void comp_drift_cur(Tile<D>& tile) = 0;

  virtual void comp_parallel_cur(Tile<D>& tile) = 0;

};


} // end of namespace ffe

