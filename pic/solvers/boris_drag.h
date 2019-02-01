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

  /// setter method to adjust drag
  //void set_drag(double drag);

  void solve(pic::Tile<D>& tile) override;

};

} // end of namespace pic
