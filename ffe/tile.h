#pragma once

#include <vector>
#include <array>

#include "../tools/mesh.h"
#include "skinny_yee.h"


namespace ffe {



/*! \brief Force-Free electrodynamics methods
 *
 */
template<std::size_t D>
class Tile : 
  virtual public fields::Tile<D>, 
  virtual public corgi::Tile<D> 
{

  public:

  // explicitly imported variables 
    
  // tile & mesh limits
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;
  using fields::Tile<D>::mesh_lengths;

  // physical parameters
  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dx;

  // main working storage grid
  using fields::Tile<D>::yee;

  // RK temporary sub-stage storages
  SkinnyYeeLattice step0;
  SkinnyYeeLattice step1;
  SkinnyYeeLattice step2;
  SkinnyYeeLattice step3;

  fields::YeeLattice& get_yee(size_t i=0) override;

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
    step0(nx,ny,nz),
    step1(nx,ny,nz),
    step2(nx,ny,nz),
    step3(nx,ny,nz)
  { }


  /// destructor
  ~Tile() override = default;

  SkinnyYeeLattice& get_step(int n);

};


} // end of namespace ffe
