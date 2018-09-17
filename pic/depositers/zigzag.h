#pragma once

#include "depositer.h"

namespace pic {

/// Linear (1st order) current depositer applying the ZigZag method
template<size_t D, size_t V>
class ZigZag :
  public virtual Depositer<D,V>
{

  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




