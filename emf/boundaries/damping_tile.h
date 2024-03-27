#pragma once

#include <type_traits>
#include <algorithm>

#include "emf/tile.h"


namespace emf {
  namespace damping {


/* \brief Damping electro-magnetic tile
 *
 */
template<
  std::size_t D,
  int S
>
class Tile : 
  virtual public emf::Tile<D>
{

  public:

  //--------------------------------------------------
  // reference field to relax tile into

  /// Electric field 
  toolbox::Mesh<float_m, 1> ex_ref;
  toolbox::Mesh<float_m, 1> ey_ref;
  toolbox::Mesh<float_m, 1> ez_ref;
  
  /// Magnetic field 
  toolbox::Mesh<float_m, 1> bx_ref;
  toolbox::Mesh<float_m, 1> by_ref;
  toolbox::Mesh<float_m, 1> bz_ref;


  /// constructor
  Tile(int nx, int ny, int nz) :
    emf::Tile<D>{nx,ny,nz},

    ex_ref{nx,ny,nz},
    ey_ref{nx,ny,nz},
    ez_ref{nx,ny,nz},

    bx_ref{nx,ny,nz},
    by_ref{nx,ny,nz},
    bz_ref{nx,ny,nz}
  { }

  //void push_e() override;
  //using Tile::push_e;
  //using Tile::push_half_b;


  //void deposit_current() override;

  using emf::Tile<D>::cfl;


  // damp field 
  void damp_fields();

  /// start index of the slope
  float_m fld1;
    
  /// end index of the slope
  float_m fld2;

  // damping parameter
  float_m ksupp = 10; 


};



  } // end on namespace damping
} // end on namespace emf
