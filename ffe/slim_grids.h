#pragma once

#include <array>
#include <vector>

#include "tools/mesh.h"
#include "definitions.h"
#include "emf/tile.h"


namespace ffe {


/// Skinny version of the full Yee lattice 
//
// Original fat version is in emf/tile.h
//
// Additionally we overload basic arithmetics to the lattice 
// to ease the implementation of RK substep combinations
class SlimGrids {

  public:

  int Nx;
  int Ny;
  int Nz;

  /// TODO: can be fiddled to be even cheaper by changing halo H=1
  // then we, however, need to implement H_i = H_j mesh copy operators

  /// Electric field 
  toolbox::Mesh<float, 0> ex;
  toolbox::Mesh<float, 0> ey;
  toolbox::Mesh<float, 0> ez;
  
  ///  Magnetic field 
  toolbox::Mesh<float, 0> bx;
  toolbox::Mesh<float, 0> by;
  toolbox::Mesh<float, 0> bz;


  // real initializer constructor
  SlimGrids(int Nx, int Ny, int Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),

    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),

    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz)
    { 
      //DEV_REGISTER
    }

  ~SlimGrids() = default;


  // basic arithmetics 
  SlimGrids& operator +=(const SlimGrids& rhs);
  SlimGrids& operator -=(const SlimGrids& rhs);
  SlimGrids& operator *=(double rhs);
  SlimGrids& operator /=(double rhs);


  void set_grids(const emf::Grids& gs);

};

SlimGrids operator +(SlimGrids lhs, const SlimGrids& rhs);
SlimGrids operator -(SlimGrids lhs, const SlimGrids& rhs);
SlimGrids operator *(SlimGrids lhs, double rhs);
SlimGrids operator /(SlimGrids lhs, double rhs);



}

//ffe::SlimGrids operator /(ffe::SlimGrids lhs, double rhs);

