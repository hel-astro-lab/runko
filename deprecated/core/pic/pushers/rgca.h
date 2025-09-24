#pragma once

#include "core/pic/pushers/pusher.h"

namespace pic {

/// Reduced guiding center approximation pusher
//
// https://arxiv.org/abs/1701.05605
//
template<size_t D, size_t V>
class rGCAPusher :
  public Pusher<D,V>
{
  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
