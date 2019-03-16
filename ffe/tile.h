#pragma once

#include <vector>
#include <array>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../em-fields/tile.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"


namespace ffe {

/// Storage for temporary lattice values
class ExtraLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// storages
  toolbox::Mesh<double, 1> ex;

  ExtraLattice(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    ex(Nx, Ny, Nz)
    { }

}




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

  ExtraLattice lattice;

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
    lattice(nx,ny,nz)
  { }


  /// destructor
  ~Tile() override = default;

  /// Compute perpendicular component of the force-free current.
  void compute_perp_curent();

  /// subtract parallel E field component to enforce E.B=0
  void subtract_parallel_e();



};



} // end of namespace ffe





