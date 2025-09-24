#pragma once

#include "core/ffe/tile.h"
#include "definitions.h"

namespace ffe {

/// Reduced FFE 4th order solver 
template<size_t D>
class rFFE4
{
  public:

  int Nx;
  int Ny;
  int Nz;

  /// interpolated arrays
  toolbox::Mesh<float, 0> bxf;
  toolbox::Mesh<float, 0> byf;
  toolbox::Mesh<float, 0> bzf;

  toolbox::Mesh<float, 0> exf;
  toolbox::Mesh<float, 0> eyf;
  toolbox::Mesh<float, 0> ezf;

  toolbox::Mesh<float, 0> rhf;


  rFFE4(int Nx, int Ny, int Nz) :
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz),
    exf(Nx, Ny, Nz),
    eyf(Nx, Ny, Nz),
    ezf(Nx, Ny, Nz),
    rhf(Nx, Ny, Nz)
  {};

  virtual ~rFFE4() = default;


  // interpolation routine
  void interpolate( 
        toolbox::Mesh<float,3>& f,
        toolbox::Mesh<float,0>& fi,
        const std::array<int,D>& in,
        const std::array<int,D>& out
      );

  /// auxiliary functions to stagger e and b
  void stagger_x_eb(emf::Grids& m);
  void stagger_y_eb(emf::Grids& m);
  void stagger_z_eb(emf::Grids& m);

  /// compute rho = div E
  void comp_rho(Tile<D>& tile);

  /// compute dE and dB
  void push_eb(Tile<D>& tile);

  /// compute and add jperp to E
  void add_jperp(Tile<D>& tile);

  /// set E_par = 0
  void remove_jpar(Tile<D>& tile);

  /// limit E < B
  void limit_e(Tile<D>& tile);


};



} // end of namespace

