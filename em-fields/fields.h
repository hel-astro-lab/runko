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

  toolbox::Mesh<Realf, 1> jx1;


  /// Kinetic energy
  toolbox::Mesh<Realf, 1> ekin;

  /// Charge density
  toolbox::Mesh<Realf, 1> rho;

  YeeLattice(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), Ny(Ny), Nz(Nz),
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

    ekin(Nx,Ny,Nz),
    rho(Nx, Ny, Nz) { }

};


/*! \brief General Plasma cell for solving Maxwell's equations
 *
 * Internally everything is stored in staggered Yee lattice.
 *
 */
class PlasmaCell : virtual public corgi::Cell {

  public:

  size_t NxMesh;
  size_t NyMesh;
  size_t NzMesh;


  // Yee lattice of plasma quantities (with 2 timesteps)
  toolbox::Rotator<YeeLattice, 1> yee;


  //--------------------------------------------------

  PlasmaCell(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh, size_t NzMesh);

  /// destructor
  ~PlasmaCell() { };


  void updateBoundaries(corgi::Node& node);

  void exchangeCurrents(corgi::Node& node);

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

  void cycleYee();

  Realf yeeDt = 1.0;
  Realf yeeDx = 1.0;

  //Realf yeeDy = 1.0;
  //Realf yeeDz = 1.0;


};



} // end of namespace fields





