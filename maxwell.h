#pragma once

#include "corgi/cell.h"
#include "corgi/corgi.h"

#include "tools/Mesh.h"

#include "dataContainer.h"



namespace maxwell {




/// Yee lattice of plasma quantities
class YeeLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// Electric field 
  toolbox::Mesh<double, 1> ex;
  toolbox::Mesh<double, 1> ey;
  toolbox::Mesh<double, 1> ez;
  
  /// Magnetic field 
  toolbox::Mesh<double, 1> bx;
  toolbox::Mesh<double, 1> by;
  toolbox::Mesh<double, 1> bz;

  /// Current vector 
  toolbox::Mesh<double, 1> jx;
  toolbox::Mesh<double, 1> jy;
  toolbox::Mesh<double, 1> jz;

  /// Charge density
  toolbox::Mesh<double, 1> rh;

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
    rh(Nx, Ny, Nz) { }

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


  // Yee lattice of plasma quantities
  datarotators::DataContainer<YeeLattice> yee;


  //--------------------------------------------------

  PlasmaCell(
      size_t i, size_t j, 
      int o,
      size_t NxG, size_t NyG,
      size_t NxMesh, size_t NyMesh);

  /// destructor
  ~PlasmaCell() { };


  void updateBoundaries(corgi::Node& node);

  void pushHalfB();
  void pushHalfB1d();
  void pushHalfB2d();
  void pushHalfB3d();

  void pushE();
  void pushE1d();
  void pushE2d();
  void pushE3d();

  void depositCurrent();


  YeeLattice& getYee();

  YeeLattice& getNewYee();


  void cycleYee();

  double dt = 0.01;
  double dx = 0.1;
  double dy = 1.0;
  double dz = 1.0;


};



} // end of namespace maxwell





