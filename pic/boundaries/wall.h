#pragma once

#include <array>
#include <vector>

#include "../../pic/tile.h"
#include "../../em-fields/damping_fields.h"


namespace pic {

template<int S>
class PicTileWall :
          public fields::PlasmaTileDamped<S>,
  virtual public pic::PicTile
{

  public:

    PicTileWall(
      size_t i, size_t j, 
      int o, 
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh
    ) : 
      corgi::Tile(i, j, o, NxG, NyG),
      fields::PlasmaTile(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
      pic::PicTile(i,j,o,NxG,NyG,NxMesh,NyMesh),
      fields::PlasmaTileDamped<S>(i,j,o,NxG, NyG, NxMesh,NyMesh, 1)
      { }

    ~PicTileWall() override = default;

    /// wall location
    using fields::PlasmaTileDamped<S>::fld1;
    using fields::PlasmaTileDamped<S>::fld2;


    // TODO: add wall movement

    // TODO: particle reflector
      
    // TODO: wall current deposition

};









}
