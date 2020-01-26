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

  int Nx;
  int Ny;
  int Nz;

  /// TODO: can be fiddled to be even cheaper by changing halo H=1
  // then we, however, need to implement H_i = H_j mesh copy operators

  /// Electric field 
  toolbox::Mesh<real_short, 3> ex;
  toolbox::Mesh<real_short, 3> ey;
  toolbox::Mesh<real_short, 3> ez;
  
  /// Magnetic field 
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
    { }

  ~SkinnyYeeLattice() = default;


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
    ex *= static_cast<real_short>(rhs);
    ey *= static_cast<real_short>(rhs);
    ez *= static_cast<real_short>(rhs);
    bx *= static_cast<real_short>(rhs);
    by *= static_cast<real_short>(rhs);
    bz *= static_cast<real_short>(rhs);

    return *this;
  }

  SkinnyYeeLattice& operator /=(double rhs) 
  {
    ex /= static_cast<real_short>(rhs);
    ey /= static_cast<real_short>(rhs);
    ez /= static_cast<real_short>(rhs);
    bx /= static_cast<real_short>(rhs);
    by /= static_cast<real_short>(rhs);
    bz /= static_cast<real_short>(rhs);

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

inline SkinnyYeeLattice operator *(SkinnyYeeLattice lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

inline SkinnyYeeLattice operator /(SkinnyYeeLattice lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}


}
