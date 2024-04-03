#pragma once

#include "core/pic/depositers/depositer.h"

namespace pic {

/// Third order current depositer applying the ZigZag method
// Uses quadratic particle shape for deposition.
template<size_t D, size_t V>
class ZigZag_3rd :
  public virtual Depositer<D,V>
{
public:
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




