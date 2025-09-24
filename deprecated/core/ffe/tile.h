#pragma once

#include <vector>
#include <array>

#include "slim_grids.h"

namespace ffe {

/*! \brief Force-Free electrodynamics methods
 *
 */
template<std::size_t D>
class Tile :
  virtual public emf::Tile<D>,
  virtual public corgi::Tile<D>
{

  public:

  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;

  using emf::Tile<D>::mesh_lengths;
  using emf::Tile<D>::grids;
  using emf::Tile<D>::cfl;


  // RK temporary sub-stage storages
  SlimGrids dF;
  SlimGrids Fn;

  /// constructor
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
     emf::Tile<D>{nx,ny,nz},
     dF{nx,ny,nz},
     Fn{nx,ny,nz}
  { }


  /// destructor
  ~Tile() override = default;


  /// update E and B
  void rk3_update(float c1, float c2, float c3);

  /// copy Y^n to Y^n-1
  void copy_eb();

};


} // end of namespace ffe
