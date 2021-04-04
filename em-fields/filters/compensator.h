#pragma once

#include "filter.h"

namespace fields {

/// 2D 3-point compensator filter from Birdsall & Langdon
template<size_t D>
class Compensator2 :
  public virtual Filter<D>
{
  public:

  using Filter<D>::Filter;

  void solve(fields::Tile<D>& tile) override;

};

} // end of namespace fields
