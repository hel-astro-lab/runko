#pragma once

#include <vector>
#include <array>

#include "skinny_yee.h"

#include"../tools/iter/allocator.h"

namespace ffe {

/*! \brief Force-Free electrodynamics methods
 *
 */
template<std::size_t D>
class Tile : 
  virtual public fields::Tile<D>, 
  virtual public corgi::Tile<D>, 
  virtual public ManagedParent
{

  public:

  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;

  using fields::Tile<D>::mesh_lengths;
  using fields::Tile<D>::yee;
  using fields::Tile<D>::cfl;


  // RK temporary sub-stage storages
  SkinnyYeeLattice dF;
  SkinnyYeeLattice Fn;

  /// constructor
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
    fields::Tile<D>{nx,ny,nz},
    dF{nx,ny,nz},
    Fn{nx,ny,nz}
  { }


  /// destructor
  ~Tile() override = default;


  /// update E and B
  void rk3_update(real_short c1, real_short c2, real_short c3);

  /// copy Y^n to Y^n-1
  void copy_eb();

};


} // end of namespace ffe
