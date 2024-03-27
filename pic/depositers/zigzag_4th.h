#pragma once

#include "pic/depositers/depositer.h"

namespace pic {

/// Second order current depositer applying the ZigZag method
// Uses triangle particle shape for deposition.
template<size_t D, size_t V>
class ZigZag_4th :
  public virtual Depositer<D,V>
{
public:
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




