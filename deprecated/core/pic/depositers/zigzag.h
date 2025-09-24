#pragma once

#include "core/pic/depositers/depositer.h"

namespace pic {

/// Linear (1st order) current depositer applying the ZigZag method
// Uses Cloud-in-cell shape (CIC)
template<size_t D, size_t V>
class ZigZag :
  public virtual Depositer<D,V>
{
public:
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




