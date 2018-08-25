#pragma once

#include <array>
#include <vector>

#include "../../pic/tile.h"
#include "../../em-fields/damping_tile.h"


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
  template< typename... Dims,
    typename = corgi::internals::enable_if_t< (sizeof...(Dims) == D) && 
               corgi::internals::are_integral<Dims...>::value, void >
  > 
  Tile(Dims... mesh_lens) :
    corgi::Tile<D>(),
    fields::Tile<D>(mesh_lens...),
    fields::damping::Tile<D,S>(mesh_lens...),
    pic::Tile<D>(mesh_lens...)
  { }

  ~Tile() override = default;

  /// wall location
  using fields::damping::Tile<D, S>::fld1;
  using fields::damping::Tile<D, S>::fld2;

  // TODO: add wall movement

  // TODO: particle reflector
      
  // TODO: wall current deposition

};







}

}
