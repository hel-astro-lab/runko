#pragma once

#include "depositer.h"

namespace pic {

/// Fourth order current depositer applying the ZigZag method
template<size_t D, size_t V>
class ZigZag_2nd :
  public virtual Depositer<D,V>
{
public:
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




