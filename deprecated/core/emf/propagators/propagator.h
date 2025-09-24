#pragma once

#include "core/emf/tile.h"
#include "definitions.h"


namespace emf {

/// General interface for Maxwell's field equation propagator
template<size_t D>
class Propagator
{
  public:

  /// time step length for more multi-stage schemes
  double dt = 1.0;

  Propagator() = default;

  virtual ~Propagator() = default;

  virtual void push_e(Tile<D>& tile) = 0;

  virtual void push_half_b(Tile<D>& tile) = 0;

};


} // end of namespace emf

