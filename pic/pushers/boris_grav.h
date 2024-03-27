#pragma once

#include "pic/pushers/pusher.h"

namespace pic {

/// Boris pusher with radiative pressure force
template<size_t D, size_t V>
class BorisPusherGrav :
  public Pusher<D,V>
{

  /// amount of drag asserted on particles
  public:
  double g0; // radiation pressure strength

  double cenx; // x location of the gravity null point 

  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
