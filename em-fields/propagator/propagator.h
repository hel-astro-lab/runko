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

  virtual void push_e(fields::Tile<D>& tile);

  virtual void push_half_b(fields::Tile<D>& tile);

};


} // end of namespace fields

