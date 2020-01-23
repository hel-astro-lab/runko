#pragma once

#include <vector>
#include <array>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../em-fields/tile.h"
#include "../tools/mesh.h"


namespace ffe {


/// Skinny version of the full Yee lattice 
//
// Original fat version is in em-fields/tile.h
//
// Additionally we overload basic arithmetics to the lattice 
// to ease the implementation of RK substep combinations
class SkinnyYeeLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// TODO: can be fiddled to be even cheaper by changing halo H=1
  // then we, however, need to implement H_i = H_j mesh copy operators

  /// Electric field 
  toolbox::Mesh<Realf, 3> ex;
  toolbox::Mesh<Realf, 3> ey;
  toolbox::Mesh<Realf, 3> ez;
  
  /// Magnetic field 
  toolbox::Mesh<Realf, 3> bx;
  toolbox::Mesh<Realf, 3> by;
  toolbox::Mesh<Realf, 3> bz;


  // real initializer constructor
  SkinnyYeeLattice(size_t Nx, size_t Ny, size_t Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),

    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),

    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz)
    { }

  virtual ~SkinnyYeeLattice() = default;


  // basic arithmetics 

  SkinnyYeeLattice& operator +=(const SkinnyYeeLattice& rhs)
  {
    ex += rhs.ex;
    ey += rhs.ey;
    ez += rhs.ez;
    bx += rhs.bx;
    by += rhs.by;
    bz += rhs.bz;

    return *this;
  }

  SkinnyYeeLattice& operator -=(const SkinnyYeeLattice& rhs) 
  {
    ex -= rhs.ex;
    ey -= rhs.ey;
    ez -= rhs.ez;
    bx -= rhs.bx;
    by -= rhs.by;
    bz -= rhs.bz;

    return *this;
  }

  SkinnyYeeLattice& operator *=(double rhs) 
  {
    ex *= rhs;
    ey *= rhs;
    ez *= rhs;
    bx *= rhs;
    by *= rhs;
    bz *= rhs;

    return *this;
  }

  SkinnyYeeLattice& operator /=(double rhs) 
  {
    ex /= rhs;
    ey /= rhs;
    ez /= rhs;
    bx /= rhs;
    by /= rhs;
    bz /= rhs;

    return *this;
  }


  // copy yee grid to skinny yee
  void set_yee(const fields::YeeLattice& yee)
  {
    ex = yee.ex;
    ey = yee.ey;
    ez = yee.ez;
    bx = yee.bx;
    by = yee.by;
    bz = yee.bz;
  }


};


inline SkinnyYeeLattice operator +(SkinnyYeeLattice lhs, const SkinnyYeeLattice& rhs)
{
  lhs += rhs;
  return lhs;
}

inline SkinnyYeeLattice operator -(SkinnyYeeLattice lhs, const SkinnyYeeLattice& rhs)
{
  lhs -= rhs;
  return lhs;
}



}
