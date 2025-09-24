#pragma once

#include <array>
#include <vector>

#include "core/pic/tile.h"
#include "core/emf/boundaries/damping_tile.h"


namespace pic {
  namespace wall {


template<size_t D, int S>
class Tile :
  virtual public emf::damping::Tile<D, S>,
  virtual public pic::Tile<D>,
  virtual public emf::Tile<D>,
  virtual public corgi::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    emf::Tile<D>(nx,ny,nz),
    emf::damping::Tile<D,S>(nx,ny,nz),
    pic::Tile<D>(nx,ny,nz)
  { }

  /// wall location
  using emf::damping::Tile<D, S>::fld1;
  using emf::damping::Tile<D, S>::fld2;

  // TODO: add wall movement

  // TODO: particle reflector
      
  // TODO: wall current deposition

};







}

}
