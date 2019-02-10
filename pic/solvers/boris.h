#pragma once

#include "pusher.h"

namespace pic {

/// Boris pusher
template<size_t D, size_t V>
class BorisPusher :
  public virtual Pusher<D,V>
{
  void push_container( pic::ParticleContainer&, double cfl);

  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic
