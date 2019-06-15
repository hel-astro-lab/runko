#pragma once

#include <vector>
#include <array>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../em-fields/tile.h"
#include "../tools/mesh.h"
#include "../tools/snapshots.h"


namespace ffe {

/// Storage for temporary lattice values
class ExtraLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// half-cell staggered B
  toolbox::Mesh<double, 1> bxf;
  toolbox::Mesh<double, 1> byf;
  toolbox::Mesh<double, 1> bzf;

  /// divergence of E
  //toolbox::Mesh<double, 1> dive;

  ExtraLattice(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz)
    //dive(Nx, Ny, Nz)
    { }

};



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
  //using fields::Tile<D>::yee;

  // explicitly import EM methods
  //using fields::Tile<D>::push_e;
  //using fields::Tile<D>::push_half_b;

  toolbox::IntegratorSnapshots<
    fields::YeeLattice> yee_snapshots;

  fields::YeeLattice& get_yee(size_t i=0) override;

  void add_yee_lattice() override;


  ExtraLattice tmp;

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
    tmp(nx,ny,nz)
  { 
  
    // preinitialize arrays for 4-stage RK
    for(size_t nstages=0; nstages < 4; nstages++)
      add_yee_lattice();

  }


  /// destructor
  ~Tile() override = default;

  /// Compute perpendicular component of the force-free current.
  void compute_perp_current();

  /// subtract parallel E field component to enforce E.B=0
  void subtract_parallel_e();


};



} // end of namespace ffe





