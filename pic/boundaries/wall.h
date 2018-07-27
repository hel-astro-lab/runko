#pragma once

#include <array>
#include <vector>

#include "../../pic/cell.h"
#include "../../em-fields/damping_fields.h"


namespace pic {

template<int S>
class PicCellWall :
          public fields::PlasmaCellDamped<S>,
  virtual public pic::PicCell
{

  public:

    PicCellWall(
      size_t i, size_t j, 
      int o, 
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh
    ) : 
      corgi::Cell(i, j, o, NxG, NyG),
      fields::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh, 1),
      pic::PicCell(i,j,o,NxG,NyG,NxMesh,NyMesh),
      fields::PlasmaCellDamped<S>(i,j,o,NxG, NyG, NxMesh,NyMesh, 1)
      { }

    ~PicCellWall() = default;

    /// wall location
    using fields::PlasmaCellDamped<S>::fld1;
    using fields::PlasmaCellDamped<S>::fld2;


    // TODO: add wall movement

    // TODO: particle reflector
      
    // TODO: wall current deposition

};









}
