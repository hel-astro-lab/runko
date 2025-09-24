#pragma once

#include "core/pic/depositers/depositer.h"

namespace pic {

/// Second order current depositer applying the Esikerpo method
template<size_t D, size_t V>
class Esikerpov_2nd :
  public virtual Depositer<D,V>
{
public:
  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic




