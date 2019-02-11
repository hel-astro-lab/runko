#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"


namespace pic {

/// General interface for PIC pushers
template<size_t D, size_t V>
class Pusher
{

  public:

  Pusher() {};

  virtual ~Pusher() = default;

  virtual void solve(pic::Tile<D>& ) = 0;

  //virtual void push_container(pic::Tile<D>& ) = 0;

  //void solve(pic::Tile<D>& tile)
  //{
  //  for(auto&& container : tile.containers)
  //    push_container(container, tile.cfl);
  //
  //}


};


} // end of namespace pic
