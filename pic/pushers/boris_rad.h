#pragma once

#include "pusher.h"

namespace pic {

/// Boris pusher with radiative pressure force
template<size_t D, size_t V>
class BorisPusherRad :
  public Pusher<D,V>
{

  /// amount of drag asserted on particles
  public:
  double drag; // radiation pressure strength

  double beam_locx; // x location of the beam front 

  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
