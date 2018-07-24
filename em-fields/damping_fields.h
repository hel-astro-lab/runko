#pragma once

#include <type_traits>
#include <algorithm>

#include "fields.h"


namespace fields {


/* \brief Damping electro-magnetic cell
 *
 */
template<int S>
class PlasmaCellDamped : 
  virtual public fields::PlasmaCell 
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

    PlasmaCellDamped(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh, size_t NzMesh);

    ~PlasmaCellDamped() { };

    //void pushE() override;
    //using PlasmaCell::pushE;
    //using PlasmaCell::pushHalfB;
      
    //void pushE2d_damped();

    void depositCurrent() override;


    using PlasmaCell::dt;
    using PlasmaCell::dx;
    using PlasmaCell::cfl;


    //-------------------------------------------------- 
    // damp field 
    //
    // Correct method is picked based on SFINAE: relies on enable_if<> returning void 

    template<int R>
    typename std::enable_if<R==-1 | R==+1>::type 
    dampFields();
 
    template<int R>
    typename std::enable_if<R==-2 | R==+2>::type 
    dampFields();
    
    // trampoline function to avoid templatizing methods
    void dampFields(){ dampFields<S>(); };


    /// start index of the slope
    Realf fld1;
      
    /// end index of the slope
    Realf fld2;

};



} // end on namespace
