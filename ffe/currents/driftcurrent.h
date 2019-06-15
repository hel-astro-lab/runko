#pragma once

#include "current.h"

namespace ffe {

/// Reduced system where only drift currents are solved
template<size_t D>
class DriftCurrent :
  public virtual Current<D>
{
  void solve(Tile<D>& tile) override;
};


} // end of namespace solve
