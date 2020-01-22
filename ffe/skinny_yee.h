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

  /// Electric field 
  toolbox::Mesh<Realf, 0> ex;
  toolbox::Mesh<Realf, 0> ey;
  toolbox::Mesh<Realf, 0> ez;
  
  /// Magnetic field 
  toolbox::Mesh<Realf, 0> bx;
  toolbox::Mesh<Realf, 0> by;
  toolbox::Mesh<Realf, 0> bz;


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
};


}
