#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"


namespace pic {

/// General interface for PIC pushers
template<size_t D, size_t V>
class Pusher
{

  public:

  Pusher() = default;

  virtual ~Pusher() = default;

  virtual void push_container(pic::ParticleContainer<D>& container, double cfl) 
  {
    // check that this is never used or that the user must know it
    assert(false);

    // initialize pointers to particle arrays
    real_prtcl* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    real_prtcl* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    for(size_t n=0; n<container.size(); n++) {
      for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n]*cfl;
    }
  }

  /// push all containers in tile
  void solve(pic::Tile<D>& tile)
  {
    for(auto&& container : tile.containers)
      push_container(container, tile.cfl);
  }


  /// push spesific containers in tile
  void solve(pic::Tile<D>& tile, int ispc)
  {
    push_container(tile.containers[ispc], tile.cfl);
  }

};


} // end of namespace pic
