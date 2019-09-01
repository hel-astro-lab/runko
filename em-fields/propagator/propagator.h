#pragma once

#include "../tile.h"
#include "../../definitions.h"


namespace fields {

/// General interface for Maxwell's field equation propagator
template<size_t D>
class Propagator
{
  public:

  Propagator() {};

  virtual ~Propagator() = default;

  virtual void push_e(Tile<D>& tile) = 0;

  virtual void push_half_b(Tile<D>& tile) = 0;

};


} // end of namespace fields

