#pragma once

#include "pusher.h"

namespace pic {

/// Boris pusher
template<size_t D, size_t V>
class BorisPusher :
  public Pusher<D,V>
{
  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
