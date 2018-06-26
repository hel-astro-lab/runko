#pragma once

#include "../corgi/cell.h"
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

  /// Current vector 
  toolbox::Mesh<Realf, 1> jx;
  toolbox::Mesh<Realf, 1> jy;
  toolbox::Mesh<Realf, 1> jz;

  /// temporary current
  toolbox::Mesh<Realf, 1> jx1;
  toolbox::Mesh<Realf, 1> jy1;
  toolbox::Mesh<Realf, 1> jz1;

  /// Charge density
  toolbox::Mesh<Realf, 1> rho;


  YeeLattice(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),

    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),

    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz),

    jx(Nx, Ny, Nz),
    jy(Nx, Ny, Nz),
    jz(Nx, Ny, Nz),

    jx1(Nx, Ny, Nz),
    jy1(Nx, Ny, Nz),
    jz1(Nx, Ny, Nz),

    rho(Nx, Ny, Nz)
    { }

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
  PlasmaMomentLattice(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),

    rho(Nx, Ny, Nz),
    mgamma(Nx, Ny, Nz),
    Vx(Nx, Ny, Nz),
    Vy(Nx, Ny, Nz),
    Vz(Nx, Ny, Nz),
    Tx(Nx, Ny, Nz),
    Ty(Nx, Ny, Nz),
    Tz(Nx, Ny, Nz),
    ekin(Nx, Ny, Nz)
  { }


};




/*! \brief General Plasma cell for solving Maxwell's equations
 *
 * Internally everything for computations are stored in 
 * staggered Yee lattice.
 *
 * Analysis/reduced data is in PlasmaMomentLattice.
 */
class PlasmaCell : virtual public corgi::Cell {

  public:

  size_t NxMesh;
  size_t NyMesh;
  size_t NzMesh;


  /// Yee lattice of plasma quantities (with 1 timestep)
  toolbox::Rotator<YeeLattice, 1> yee;

  /// species specific analysis results
  std::vector<PlasmaMomentLattice> analysis;


  //--------------------------------------------------

  PlasmaCell(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh, size_t NzMesh);

  /// destructor
  ~PlasmaCell() { };


  void updateBoundaries(corgi::Node& node);
  void updateBoundaries2D(corgi::Node& node);

  void exchangeCurrents(corgi::Node& node);
  void exchangeCurrents2D(corgi::Node& node);

  virtual void pushHalfB();
  void pushHalfB1d();
  void pushHalfB2d();
  void pushHalfB3d();

  virtual void pushE();
  void pushE1d();
  void pushE2d();
  void pushE3d();

  virtual void depositCurrent();

  YeeLattice& getYee(size_t i=0);

  PlasmaMomentLattice& getAnalysis(size_t i);

  void cycleYee();

  void cycleCurrent();

  Realf yeeDt = 1.0;
  Realf yeeDx = 1.0;

  //Realf yeeDy = 1.0;
  //Realf yeeDz = 1.0;

  void addAnalysisSpecies();
};



} // end of namespace fields





