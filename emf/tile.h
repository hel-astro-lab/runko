#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../tools/mesh.h"
#include "../definitions.h"
#include"../tools/iter/allocator.h"

namespace emf {
  namespace mpi = mpi4cpp::mpi;

/// Yee lattice of plasma quantities
class YeeLattice
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
  YeeLattice()  = default;

  // real initializer constructor
  YeeLattice(int Nx, int Ny, int Nz) : 
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
  YeeLattice(YeeLattice& other) = default;

  YeeLattice(const YeeLattice& other) = default;

  // move constructor
  YeeLattice(YeeLattice&& other) noexcept :
      //YeeLattice() // initialize via default constructor, C++11 only
      YeeLattice{other.Nx, other.Ny, other.Nz} // initialize via allocating constructor
  {
    swap(*this, other);
  }


  // public swap for efficient memory management
  friend void swap(YeeLattice& first, YeeLattice& second)
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
  YeeLattice& operator=(YeeLattice other) 
  {
    swap(*this, other);
    return *this;
  }

  ~YeeLattice()
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
  YeeLattice yee;

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
     yee{nx, ny, nz}
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
  virtual YeeLattice& get_yee(int /*i*/ =0){ return this->yee; }

  /// Get current time snapshot of Yee lattice
  virtual YeeLattice& get_yee2(){ return this->yee; }

  /// Set current time snapshot of Yee lattice
  virtual void set_yee(YeeLattice& v){ this->yee = v; }

  // not needed anymore; also undefined behavior since this can be deleted
  //virtual std::shared_ptr<YeeLattice> get_yeeptr(){ return std::shared_ptr<YeeLattice>(&yee); }

  //FIXME: naming is wrong because of pybind: it can not deduce correct function
  //       with only constness difference.
  virtual const YeeLattice& get_const_yee(int i=0) const { return this->yee; }

  virtual void clear_current();

  std::vector<mpi::request> 
  send_data( mpi::communicator& /*comm*/, int dest, int mode, int tag) override;

  std::vector<mpi::request> 
  recv_data( mpi::communicator& /*comm*/, int orig, int mode, int tag) override;

};


} // end of namespace emf





