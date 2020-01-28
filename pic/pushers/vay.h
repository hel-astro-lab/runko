#pragma once

#include "pusher.h"

namespace pic {

/// Vay pusher
template<size_t D, size_t V>
class VayPusher :
  public Pusher<D,V>
{
  void push_container( pic::ParticleContainer<D>&, double cfl) override;
};

} // end of namespace pic
