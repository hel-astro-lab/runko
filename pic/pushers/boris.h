#pragma once

#include "pusher.h"

namespace pic {

/// Boris pusher
template<size_t D, size_t V>
class BorisPusher :
  public Pusher<D,V>
{
  public: // again why is this not public ?
  void push_container( pic::ParticleContainer<D>& /*container*/, double cfl) override;
};

} // end of namespace pic
