#pragma once

#include <array>
#include <vector>

#include "../tools/mesh.h"
#include "../definitions.h"
#include "../em-fields/tile.h"


namespace ffe {


/// Skinny version of the full Yee lattice 
//
// Original fat version is in em-fields/tile.h
//
// Additionally we overload basic arithmetics to the lattice 
// to ease the implementation of RK substep combinations
class SkinnyYeeLattice {

  public:

  int Nx;
  int Ny;
  int Nz;

  /// TODO: can be fiddled to be even cheaper by changing halo H=1
  // then we, however, need to implement H_i = H_j mesh copy operators

  /// Electric field 
  toolbox::Mesh<real_short, 3> ex;
  toolbox::Mesh<real_short, 3> ey;
  toolbox::Mesh<real_short, 3> ez;
  
  ///  Magnetic field 
  toolbox::Mesh<real_short, 3> bx;
  toolbox::Mesh<real_short, 3> by;
  toolbox::Mesh<real_short, 3> bz;


  // real initializer constructor
  SkinnyYeeLattice(int Nx, int Ny, int Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),

    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),

    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz)
    { 
      DEV_REGISTER
    }

  ~SkinnyYeeLattice() = default;


  // basic arithmetics 
  SkinnyYeeLattice& operator +=(const SkinnyYeeLattice& rhs);
  SkinnyYeeLattice& operator -=(const SkinnyYeeLattice& rhs);
  SkinnyYeeLattice& operator *=(double rhs);
  SkinnyYeeLattice& operator /=(double rhs);


  void set_yee(const fields::YeeLattice& yee);

};

SkinnyYeeLattice operator +(SkinnyYeeLattice lhs, const SkinnyYeeLattice& rhs);
SkinnyYeeLattice operator -(SkinnyYeeLattice lhs, const SkinnyYeeLattice& rhs);
SkinnyYeeLattice operator *(SkinnyYeeLattice lhs, double rhs);
SkinnyYeeLattice operator /(SkinnyYeeLattice lhs, double rhs);



}

//ffe::SkinnyYeeLattice operator /(ffe::SkinnyYeeLattice lhs, double rhs);

