#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "corgi/tile.h"
#include "corgi/corgi.h"
#include "tools/mesh.h"
#include "definitions.h"
#include "tools/iter/allocator.h"

namespace emf {
  namespace mpi = mpi4cpp::mpi;

/// Yee lattice of plasma quantities
class Grids
{

  public:

  int Nx;
  int Ny;
  int Nz;

  /// Electric field 
  toolbox::Mesh<float_m, 3> ex;
  toolbox::Mesh<float_m, 3> ey;
  toolbox::Mesh<float_m, 3> ez;
  
  /// Magnetic field 
  toolbox::Mesh<float_m, 3> bx;
  toolbox::Mesh<float_m, 3> by;
  toolbox::Mesh<float_m, 3> bz;
    
  /// Charge density
  toolbox::Mesh<float_m, 3> rho;

  /// Current vector 
  toolbox::Mesh<float_m, 3> jx;
  toolbox::Mesh<float_m, 3> jy;
  toolbox::Mesh<float_m, 3> jz;



  // default empty constructor
  Grids()  = default;

  // real initializer constructor
  Grids(int Nx, int Ny, int Nz) : 
    Nx{Nx}, Ny{Ny}, Nz{Nz},
    ex{Nx, Ny, Nz},
    ey{Nx, Ny, Nz},
    ez{Nx, Ny, Nz},
    bx{Nx, Ny, Nz},
    by{Nx, Ny, Nz},
    bz{Nx, Ny, Nz},
    rho{Nx, Ny, Nz},
    jx{Nx, Ny, Nz},
    jy{Nx, Ny, Nz},
    jz{Nx, Ny, Nz}
  { 
    //DEV_REGISTER
  }

  // copy ctor
  Grids(Grids& other) = default;

  Grids(const Grids& other) = default;

  // move constructor
  Grids(Grids&& other) noexcept :
      //Grids() // initialize via default constructor, C++11 only
      Grids{other.Nx, other.Ny, other.Nz} // initialize via allocating constructor
  {
    swap(*this, other);
  }


  // public swap for efficient memory management
  friend void swap(Grids& first, Grids& second)
  {
    using std::swap;
    swap(first.Nx,  second.Nx);
    swap(first.Ny,  second.Ny);
    swap(first.Nz,  second.Nz);
    swap(first.ex , second.ex);
    swap(first.ey , second.ey);
    swap(first.ez , second.ez);
    swap(first.bx , second.bx);
    swap(first.by , second.by);
    swap(first.bz , second.bz);
    swap(first.rho, second.rho);
    swap(first.jx , second.jx);
    swap(first.jy , second.jy);
    swap(first.jz , second.jz);


  }

  // copy-and-swap algorithm
  Grids& operator=(Grids other) 
  {
    swap(*this, other);
    return *this;
  }

  ~Grids()
  {

  };

};



/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in 
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : 
  virtual public corgi::Tile<D>,
  virtual public ManagedParent
{

  public:

  /// size of grid inside tile
  std::array<int, 3> mesh_lengths;

  /// Yee lattice of plasma quantities (with 1 timestep)
  Grids grids;

  /// explicitly show that we import tile limits from base class 
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;

  /// CFL number (corresponds to simulation light speed c)
  double cfl;

  //--------------------------------------------------
  // constructor with internal mesh dimensions
  Tile(int nx, int ny, int nz) :
     corgi::Tile<D>(),
     mesh_lengths {{nx, ny, nz}},
     grids{nx, ny, nz}
  {
    if (D == 1) assert(ny == 1 && nz == 1);
    if (D == 2) assert(nz == 1);

  }

  ~Tile() override = default;

  //--------------------------------------------------

  virtual void update_boundaries(  
          corgi::Grid<D>& grid,
          std::vector<int> iarr = {0,1,2}
          );

  virtual void exchange_currents(  corgi::Grid<D>& grid);

  virtual void deposit_current();

  /// Get current time snapshot of Yee lattice
  virtual Grids& get_grids(int /*i*/ =0){ return this->grids; }

  /// Get current time snapshot of Yee lattice
  virtual Grids& get_grids2(){ return this->grids; }

  /// Set current time snapshot of Yee lattice
  virtual void set_grids(Grids& v){ this->grids = v; }

  // not needed anymore; also undefined behavior since this can be deleted
  //virtual std::shared_ptr<Grids> get_grids_ptr(){ return std::shared_ptr<Grids>(&grids); }

  //FIXME: naming is wrong because of pybind: it can not deduce correct function
  //       with only constness difference.
  virtual const Grids& get_const_grids(int /*i*/=0) const { return this->grids; }

  virtual void clear_current();

  std::vector<mpi::request> 
  send_data( mpi::communicator& /*comm*/, int dest, int mode, int tag) override;

  std::vector<mpi::request> 
  recv_data( mpi::communicator& /*comm*/, int orig, int mode, int tag) override;

};


} // end of namespace emf





