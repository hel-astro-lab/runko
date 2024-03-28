#pragma once

#include "core/pic/pushers/pusher.h"

namespace pic {

/// Higuera-Cary pusher
//
// https://arxiv.org/pdf/1701.05605.pdf
//
template<size_t D, size_t V>
class HigueraCaryPusher :
  public Pusher<D,V>
{
  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
