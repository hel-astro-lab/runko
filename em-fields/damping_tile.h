#pragma once

#include <type_traits>
#include <algorithm>

#include "tile.h"


namespace fields {
  namespace damping {


/* \brief Damping electro-magnetic tile
 *
 */
template<
  std::size_t D,
  int S
>
class Tile : 
  virtual public fields::Tile<D>
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


  /// constructor
  template< typename... Dims,
    typename = corgi::internals::enable_if_t< (sizeof...(Dims) == D) && 
               corgi::internals::are_integral<Dims...>::value, void >
  > 
  Tile(Dims... mesh_lens) :
     corgi::Tile<D>(),
    fields::Tile<D>(mesh_lens...),

    ex_ref(mesh_lens...),
    ey_ref(mesh_lens...),
    ez_ref(mesh_lens...),

    bx_ref(mesh_lens...),
    by_ref(mesh_lens...),
    bz_ref(mesh_lens...)
  { }


  ~Tile() override = default;

  //void pushE() override;
  //using Tile::pushE;
  //using Tile::pushHalfB;
    
  //void pushE2d_damped();

  void depositCurrent() override;

  //using fields::Tile<D>::dt;
  using fields::Tile<D>::dx;
  using fields::Tile<D>::cfl;

  // damp field 
  void dampFields();

  /// start index of the slope
  Realf fld1;
    
  /// end index of the slope
  Realf fld2;

};



  } // end on namespace damping
} // end on namespace fields
