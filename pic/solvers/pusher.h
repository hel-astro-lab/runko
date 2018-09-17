#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"


namespace pic {

/// General interface for PIC pushers
template<size_t D, size_t V>
class Pusher
{

  public:

  Pusher() {};

  virtual ~Pusher() = default;

  virtual void solve(pic::Tile<D>& ) = 0;

};


} // end of namespace pic
