#pragma once

#include "DenseGrid.h"



namespace maxwell {

/*! \brief General Plasma cell for solving Maxwell's equations
 *
 * Internally everything is stored in staggered Yee lattice.
 *
 */
class PlasmaCell {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;


  // Yee lattice of plasma quantities
  //--------------------------------------------------
  
  /// Electric field 
  toolbox::DenseGrid<double> ex;
  toolbox::DenseGrid<double> ey;
  toolbox::DenseGrid<double> ez;
  
  /// Magnetic field 
  toolbox::DenseGrid<double> bx;
  toolbox::DenseGrid<double> by;
  toolbox::DenseGrid<double> bz;

  /// Current vector 
  toolbox::DenseGrid<double> jx;
  toolbox::DenseGrid<double> jy;
  toolbox::DenseGrid<double> jz;

  /// Charge density
  toolbox::DenseGrid<double> rho;

  //--------------------------------------------------

  PlasmaCell(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), Ny(Ny), Nz(Nz),
    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),
    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz),
    jx(Nx, Ny, Nz),
    jy(Nx, Ny, Nz),
    jz(Nx, Ny, Nz),
    rho(Nx, Ny, Nz) { }


  /// Update B field with a half step
  void pushHalfB();

  /// Update E field with full step
  void pushE();


};



} // end of namespace maxwell





