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

  // external emf that are added to the pusher Lorentz force
  double bx_ext = 0.0;
  double by_ext = 0.0;
  double bz_ext = 0.0;

  double ex_ext = 0.0;
  double ey_ext = 0.0;
  double ez_ext = 0.0;

  // external field getter functions; can be overridden in derived pusher classes to depend on position
  DEVCALLABLE virtual double get_bx_ext(double /*x*/, double /*y*/, double /*z*/) { return bx_ext; };
  DEVCALLABLE virtual double get_by_ext(double /*x*/, double /*y*/, double /*z*/) { return by_ext; };
  DEVCALLABLE virtual double get_bz_ext(double /*x*/, double /*y*/, double /*z*/) { return bz_ext; };

  DEVCALLABLE virtual double get_ex_ext(double /*x*/, double /*y*/, double /*z*/) { return ex_ext; };
  DEVCALLABLE virtual double get_ey_ext(double /*x*/, double /*y*/, double /*z*/) { return ey_ext; };
  DEVCALLABLE virtual double get_ez_ext(double /*x*/, double /*y*/, double /*z*/) { return ez_ext; };


  virtual void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& /*tile*/
          ) 
  {
    // check that this is never used or that the user must know it
    assert(false);

    // initialize pointers to particle arrays
    float_p* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    float_p* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    for(size_t n=0; n<container.size(); n++) {
      for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n];
    }
  }

  /// push all containers in tile
  void solve(pic::Tile<D>& tile)
  {
    for(auto&& container : tile.containers)
      push_container(container, tile);
  }


  /// push spesific containers in tile
  void solve(pic::Tile<D>& tile, int ispc)
  {
    push_container(tile.containers[ispc], tile);
  }

};


} // end of namespace pic
