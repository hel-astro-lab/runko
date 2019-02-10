#pragma once

#include "pusher.h"

namespace pic {

/// Boris pusher with drag force
template<size_t D, size_t V>
class BorisPusherDrag :
  public virtual Pusher<D,V>
{

  /// amount of drag asserted on particles
  public:
  double drag;

  void push_container( pic::ParticleContainer&, double cfl);

  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic
