#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../definitions.h"


namespace fields {
  namespace mpi = mpi4cpp::mpi;


/// Yee lattice of plasma quantities
class YeeLattice
// : public std::enable_shared_from_this<YeeLattice>
{

  public:

  int Nx;
  int Ny;
  int Nz;

  /// Electric field 
  toolbox::Mesh<real_short, 3> ex;
  toolbox::Mesh<real_short, 3> ey;
  toolbox::Mesh<real_short, 3> ez;
  
  /// Magnetic field 
  toolbox::Mesh<real_short, 3> bx;
  toolbox::Mesh<real_short, 3> by;
  toolbox::Mesh<real_short, 3> bz;
    
  /// Charge density
  toolbox::Mesh<real_short, 1> rho;

  /// Current vector 
  toolbox::Mesh<real_short, 3> jx;
  toolbox::Mesh<real_short, 3> jy;
  toolbox::Mesh<real_short, 3> jz;

  // default empty constructor
  YeeLattice()  = default;

  // real initializer constructor
  YeeLattice(int Nx, int Ny, int Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),
    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz),
    rho(Nx, Ny, Nz),
    jx(Nx, Ny, Nz),
    jy(Nx, Ny, Nz),
    jz(Nx, Ny, Nz)
  {
    Nx = Nx; 
    Ny = Ny; 
    Nz = Nz; 
  }

  // copy ctor
  YeeLattice(YeeLattice& other) :
    Nx(other.Nx),
    Ny(other.Ny),
    Nz(other.Nz),
    ex(other.ex),
    ey(other.ey),
    ez(other.ez),
    bx(other.bx),
    by(other.by),
    bz(other.bz),
    rho(other.rho),
    jx(other.jx),
    jy(other.jy),
    jz(other.jz)
  {
    Nx = other.Nx; 
    Ny = other.Ny; 
    Nz = other.Nz; 
  }

  YeeLattice(const YeeLattice& other) :
    Nx(other.Nx),
    Ny(other.Ny),
    Nz(other.Nz),
    ex(other.ex),
    ey(other.ey),
    ez(other.ez),
    bx(other.bx),
    by(other.by),
    bz(other.bz),
    rho(other.rho),
    jx(other.jx),
    jy(other.jy),
    jz(other.jz)
  { 
    Nx = other.Nx; 
    Ny = other.Ny; 
    Nz = other.Nz; 
  }

  // move constructor
  YeeLattice(YeeLattice&& other)
      : YeeLattice() // initialize via default constructor, C++11 only
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

  ~YeeLattice() = default;

};


/// Lattice to hold plasma moment values of separate species
class PlasmaMomentLattice 
//  : public std::enable_shared_from_this<PlasmaMomentLattice>
{

  public:

  int Nx;
  int Ny;
  int Nz;

  // density
  toolbox::Mesh<real_short, 0> rho;

  /// energy density
  toolbox::Mesh<real_short, 0> edens;

  /// temperature
  toolbox::Mesh<real_short, 0> temp;

  /// (mean) flow speed
  toolbox::Mesh<real_short, 0> Vx;
  toolbox::Mesh<real_short, 0> Vy;
  toolbox::Mesh<real_short, 0> Vz;

  /// momentum density
  toolbox::Mesh<real_short, 0> momx;
  toolbox::Mesh<real_short, 0> momy;
  toolbox::Mesh<real_short, 0> momz;

  /// pressure
  toolbox::Mesh<real_short, 0> pressx;
  toolbox::Mesh<real_short, 0> pressy;
  toolbox::Mesh<real_short, 0> pressz;

  /// shear
  toolbox::Mesh<real_short, 0> shearxy;
  toolbox::Mesh<real_short, 0> shearxz;
  toolbox::Mesh<real_short, 0> shearyz;

  // default empty constructor
  PlasmaMomentLattice() {};


  /// constructor; initializes Mesh objects at the same time
  PlasmaMomentLattice(int Nx_in, int Ny_in, int Nz_in) : 
    Nx(Nx_in), Ny(Ny_in), Nz(Nz_in),

    rho(Nx, Ny, Nz),
    edens(Nx, Ny, Nz),
    temp(Nx, Ny, Nz),
    Vx( Nx, Ny, Nz ),
    Vy( Nx, Ny, Nz ),
    Vz( Nx, Ny, Nz ),
    momx( Nx, Ny, Nz ),
    momy( Nx, Ny, Nz ),
    momz( Nx, Ny, Nz ),
    pressx( Nx, Ny, Nz ),
    pressy( Nx, Ny, Nz ),
    pressz( Nx, Ny, Nz ),
    shearxy( Nx, Ny, Nz ),
    shearxz( Nx, Ny, Nz ),
    shearyz( Nx, Ny, Nz )
  { }

  //virtual ~PlasmaMomentLattice() = default;
  ~PlasmaMomentLattice() = default;

  // clear all internal storages
  void clear() 
  {
    rho    .clear();
    edens  .clear();
    temp   .clear();
    Vx     .clear();   
    Vy     .clear();  
    Vz     .clear();     
    momx   .clear();      
    momy   .clear();       
    momz   .clear();       
    pressx .clear();        
    pressy .clear();       
    pressz .clear();       
    shearxy.clear();
    shearxz.clear();
    shearyz.clear();  
  }
};




/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in 
 * staggered Yee lattice.
 *
 * Analysis/reduced data is in PlasmaMomentLattice.
 */
template<std::size_t D>
class Tile : 
  virtual public corgi::Tile<D>
{

  public:

  /// size of grid inside tile
  std::array<int, 3> mesh_lengths;

  /// Yee lattice of plasma quantities (with 1 timestep)
  //toolbox::Rotator<YeeLattice, 1> yee;
  //YeeLattice yee;
  std::vector<YeeLattice> yee;

  /// species-specific analysis results
  std::vector<PlasmaMomentLattice> analysis;

  /// explicitly show that we import tile limits from base class 
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;


  /// CFL number (corresponds to simulation light speed c)
  real_long cfl;

  /// grid size (assuming cubical cells)
  //real_long dx = 1.0;

  //--------------------------------------------------
  // constructor with internal mesh dimensions
  Tile(int nx, int ny, int nz) :
    mesh_lengths {{nx, ny, nz}}
    //yee(nx, ny, nz)
  {
    if (D == 1) assert(ny == 1 && nz == 1);
    if (D == 2) assert(nz == 1);

    // initialize one Yee lattice into the grid 
    // TODO: into data rotator?

    add_yee_lattice();
  }

  // avoid copies; TODO: is this needed?
  //Tile(Tile& ) = delete;

  //~Tile() = default;
  virtual ~Tile() = default;

  //--------------------------------------------------

  virtual void update_boundaries(  corgi::Grid<D>& grid);

  virtual void exchange_currents(  corgi::Grid<D>& grid);

  virtual void deposit_current();

  virtual YeeLattice& get_yee(int i=0);

  virtual YeeLattice& get_yee2();

  virtual void set_yee(YeeLattice& v);

  virtual std::shared_ptr<YeeLattice> get_yeeptr();

  virtual PlasmaMomentLattice& get_analysis(int i);

  //FIXME: naming is wrong because of pybind: it can not deduce correct function
  //       with only constness difference.
  virtual const YeeLattice& get_const_yee(int i=0) const; 

  virtual const PlasmaMomentLattice& get_const_analysis(int i) const;

  void cycle_yee();

  void cycle_current();

  void clear_current();

  void add_yee_lattice();

  void add_analysis_species();

  virtual std::vector<mpi::request> 
  send_data( mpi::communicator&, int orig, int mode, int tag) override;

  virtual std::vector<mpi::request> 
  recv_data( mpi::communicator&, int dest, int mode, int tag) override;
};



} // end of namespace fields





