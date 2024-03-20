#pragma once

#include <array>
#include <vector>

#include "../../pic/tile.h"
#include "../../emf/boundaries/damping_tile.h"


namespace pic {
  namespace wall {


template<size_t D, int S>
class Tile :
  virtual public fields::damping::Tile<D, S>,
  virtual public pic::Tile<D>,
  virtual public fields::Tile<D>,
  virtual public corgi::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
    fields::damping::Tile<D,S>(nx,ny,nz),
    pic::Tile<D>(nx,ny,nz)
  { }

  /// wall location
  using fields::damping::Tile<D, S>::fld1;
  using fields::damping::Tile<D, S>::fld2;

  // TODO: add wall movement

  // TODO: particle reflector
      
  // TODO: wall current deposition

};







}

}
