#pragma once

#include <vector>
#include <array>
#include <mpi4cpp/mpi.h>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../em-fields/tile.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../definitions.h"


namespace ffe {
  namespace mpi = mpi4cpp::mpi;


/*! \brief Force-Free electrodynamics methods
 *
 */
template<std::size_t D>
class Tile : 
  virtual public fields::Tile<D>, 
  virtual public corgi::Tile<D> 
{

  public:

  /// explicitly variables import what we need
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;

  using fields::Tile<D>::mesh_lengths;
  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dx;
  using fields::Tile<D>::yee;


  // explicitly import EM methods
  using fields::Tile<D>::push_e();
  using fields::Tile<D>::push_half_b();


  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz)
  { }


  /// destructor
  ~Tile() override = default;


  /// Compute perpendicular component of the force-free current.
  void compute_perp_curent();




};



} // end of namespace ffe





