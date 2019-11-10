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
class YeeLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// Electric field 
  toolbox::Mesh<Realf, 3> ex;
  toolbox::Mesh<Realf, 3> ey;
  toolbox::Mesh<Realf, 3> ez;
  
  /// Magnetic field 
  toolbox::Mesh<Realf, 3> bx;
  toolbox::Mesh<Realf, 3> by;
  toolbox::Mesh<Realf, 3> bz;
    
  /// Charge density
  toolbox::Mesh<Realf, 1> rho;

  /// Current vector 
  toolbox::Mesh<Realf, 3> jx;
  toolbox::Mesh<Realf, 3> jy;
  toolbox::Mesh<Realf, 3> jz;

  /// temporary current
  toolbox::Mesh<Realf, 3> jx1;
  toolbox::Mesh<Realf, 3> jy1;
  toolbox::Mesh<Realf, 3> jz1;


  // default empty constructor
  //YeeLattice() {};

  // real initializer constructor
  YeeLattice(size_t Nx, size_t Ny, size_t Nz) : 
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
    jz(Nx, Ny, Nz),

    jx1(Nx, Ny, Nz),
    jy1(Nx, Ny, Nz),
    jz1(Nx, Ny, Nz)
    { }

  virtual ~YeeLattice() = default;
};


/// Lattice to hold plasma moment values of separate species
class PlasmaMomentLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  // density
  toolbox::Mesh<Realf, 0> rho;

  /// energy density
  toolbox::Mesh<Realf, 0> edens;

  /// temperature
  toolbox::Mesh<Realf, 0> temp;

  /// (mean) flow speed
  toolbox::Mesh<Realf, 0> Vx;
  toolbox::Mesh<Realf, 0> Vy;
  toolbox::Mesh<Realf, 0> Vz;

  /// momentum density
  toolbox::Mesh<Realf, 0> momx;
  toolbox::Mesh<Realf, 0> momy;
  toolbox::Mesh<Realf, 0> momz;

  /// pressure
  toolbox::Mesh<Realf, 0> pressx;
  toolbox::Mesh<Realf, 0> pressy;
  toolbox::Mesh<Realf, 0> pressz;

  /// shear
  toolbox::Mesh<Realf, 0> shearxy;
  toolbox::Mesh<Realf, 0> shearxz;
  toolbox::Mesh<Realf, 0> shearyz;


  /// constructor; initializes Mesh objects at the same time
  PlasmaMomentLattice(size_t Nx_in, size_t Ny_in, size_t Nz_in) : 
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

  virtual ~PlasmaMomentLattice() = default;

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
  std::array<size_t, 3> mesh_lengths;

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
  Realf cfl;

  /// grid size (assuming cubical cells)
  Realf dx = 1.0;

  //--------------------------------------------------
  // constructor with internal mesh dimensions
  Tile(size_t nx, size_t ny, size_t nz) :
    mesh_lengths {{nx, ny, nz}}
  {
    if (D == 1) assert(ny == 1 && nz == 1);
    if (D == 2) assert(nz == 1);

    // initialize one Yee lattice into the grid 
    // TODO: into data rotator?
    add_yee_lattice();
  }

  // avoid copies; TODO: is this needed?
  Tile(Tile& ) = delete;

  ~Tile() = default;

  //--------------------------------------------------

  virtual void update_boundaries(  corgi::Grid<D>& grid);

  virtual void exchange_currents(  corgi::Grid<D>& grid);

  virtual void deposit_current();

  virtual YeeLattice& get_yee(size_t i=0);

  virtual PlasmaMomentLattice& get_analysis(size_t i);

  //FIXME: naming is wrong because of pybind: it can not deduce correct function
  //       with only constness difference.
  virtual const YeeLattice& get_const_yee(size_t i=0) const; 

  virtual const PlasmaMomentLattice& get_const_analysis(size_t i) const;

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





