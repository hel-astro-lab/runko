#pragma once

#include "core/pic/pushers/pusher.h"

namespace pic {

/// Simple photon pusher
template<size_t D, size_t V>
class PhotonPusher :
  public Pusher<D,V>
{
  public: // again why is this not public ?
  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;

};

} // end of namespace pic
