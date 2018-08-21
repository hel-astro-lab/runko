#pragma once

#include <type_traits>
#include <algorithm>

#include "fields.h"


namespace fields {


/* \brief Damping electro-magnetic tile
 *
 */
template<int S>
class PlasmaTileDamped : 
  virtual public fields::PlasmaTile 
{

  public:

    //--------------------------------------------------
    // reference field to relax tile into

    /// Electric field 
    toolbox::Mesh<Realf, 1> ex_ref;
    toolbox::Mesh<Realf, 1> ey_ref;
    toolbox::Mesh<Realf, 1> ez_ref;
    
    /// Magnetic field 
    toolbox::Mesh<Realf, 1> bx_ref;
    toolbox::Mesh<Realf, 1> by_ref;
    toolbox::Mesh<Realf, 1> bz_ref;

    PlasmaTileDamped(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh, size_t NzMesh);

    ~PlasmaTileDamped() override = default;

    //void pushE() override;
    //using PlasmaTile::pushE;
    //using PlasmaTile::pushHalfB;
      
    //void pushE2d_damped();

    void depositCurrent() override;


    using PlasmaTile::dt;
    using PlasmaTile::dx;
    using PlasmaTile::cfl;


    // damp field 
    void dampFields();

    /// start index of the slope
    Realf fld1;
      
    /// end index of the slope
    Realf fld2;

};



} // end on namespace
