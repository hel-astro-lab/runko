#pragma once

#include "../corgi/tile.h"
#include "../corgi/corgi.h"

#include "../tools/mesh.h"
#include "../tools/rotator.h"

#include "../definitions.h"


namespace fields {


/// Yee lattice of plasma quantities
class YeeLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// Electric field 
  toolbox::Mesh<Realf, 1> ex;
  toolbox::Mesh<Realf, 1> ey;
  toolbox::Mesh<Realf, 1> ez;
  
  /// Magnetic field 
  toolbox::Mesh<Realf, 1> bx;
  toolbox::Mesh<Realf, 1> by;
  toolbox::Mesh<Realf, 1> bz;
    
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

  ~YeeLattice() = default;

};


/// Lattice to hold plasma moment values of separate species
class PlasmaMomentLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  // density
  toolbox::Mesh<Realf, 0> rho;

  /// mean gamma
  toolbox::Mesh<Realf, 0> mgamma;

  /// (mean) flow speed
  toolbox::Mesh<Realf, 0> Vx;
  toolbox::Mesh<Realf, 0> Vy;
  toolbox::Mesh<Realf, 0> Vz;

  /// Non-relativistic temperature
  toolbox::Mesh<Realf, 0> Tx;
  toolbox::Mesh<Realf, 0> Ty;
  toolbox::Mesh<Realf, 0> Tz;

  /// Kinetic energy
  toolbox::Mesh<Realf, 0> ekin;


  /// constructor; initializes Mesh objects at the same time
  PlasmaMomentLattice(size_t Nx_in, size_t Ny_in, size_t Nz_in) : 
    Nx(Nx_in), Ny(Ny_in), Nz(Nz_in),
    rho(Nx, Ny, Nz),
    mgamma(Nx, Ny, Nz),
    Vx( Nx, Ny, Nz ),
    Vy( Nx, Ny, Nz ),
    Vz( Nx, Ny, Nz ),
    Tx( Nx, Ny, Nz ),
    Ty( Nx, Ny, Nz ),
    Tz( Nx, Ny, Nz ),
    ekin(Nx, Ny, Nz)
  { }


  ~PlasmaMomentLattice() = default;

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
  toolbox::Rotator<YeeLattice, 1> yee;

  /// species specific analysis results
  std::vector<PlasmaMomentLattice> analysis;

  /// explicitly show that we import tile limits from base class 
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;


  /// CFL number (corresponds to simulation light speed c)
  Realf cfl;

  /// grid size (assuming cubical cells)
  Realf dx = 1.0;

  /// simulation timestep
  // TODO: removable?
  //Realf dt = 1.0;


  //--------------------------------------------------
  // constructor with internal mesh dimensions
  Tile(size_t nx, size_t ny, size_t nz) :
    mesh_lengths {{nx, ny, nz}}
  {
    // initialize one Yee lattice into the grid (into data rotator)
    addYeeLattice();
  }


  /// destructor
  virtual ~Tile() override = default;


  Tile(Tile& ) = delete;


  //--------------------------------------------------

  virtual void updateBoundaries(  corgi::Node<D>& node);

  virtual void exchangeCurrents(  corgi::Node<D>& node);

  virtual void pushHalfB();

  virtual void pushE();

  virtual void depositCurrent();

  virtual YeeLattice& getYee(size_t i=0);

  virtual PlasmaMomentLattice& getAnalysis(size_t i);

  void cycleYee();

  void cycleCurrent();

  void addYeeLattice();

  void addAnalysisSpecies();
};



} // end of namespace fields





