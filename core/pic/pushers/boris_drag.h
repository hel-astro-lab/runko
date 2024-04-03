#pragma once

#include "core/pic/pushers/pusher.h"

namespace pic {

/// Boris pusher with drag force
template<size_t D, size_t V>
class BorisPusherDrag :
  public Pusher<D,V>
{

  /// amount of drag asserted on particles
  public:
  double drag; // = gamma_rad^-2
  double temp; // = 3\Theta

  double freezing_factor = 1.0; // [0,1] parameter to define how much particles move

  // Klein-Nishina cross-section
  double kn(double x);

  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
