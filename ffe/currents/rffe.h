#pragma once

#include "../tile.h"
#include "../../definitions.h"

namespace ffe {

/// Combined RK3 reduced FFE 2nd order solver 
template<size_t D>
class rFFE2
{
  public:

  int Nx;
  int Ny;
  int Nz;

  /// interpolated arrays
  toolbox::Mesh<real_short, 0> bxf;
  toolbox::Mesh<real_short, 0> byf;
  toolbox::Mesh<real_short, 0> bzf;

  toolbox::Mesh<real_short, 0> exf;
  toolbox::Mesh<real_short, 0> eyf;
  toolbox::Mesh<real_short, 0> ezf;

  toolbox::Mesh<real_short, 0> rhf;


  rFFE2(int Nx, int Ny, int Nz) :
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz),
    exf(Nx, Ny, Nz),
    eyf(Nx, Ny, Nz),
    ezf(Nx, Ny, Nz),
    rhf(Nx, Ny, Nz)
  {
    DEV_REGISTER
  };

  virtual ~rFFE2() = default;


  // interpolation routine
  void interpolate( 
        toolbox::Mesh<real_short,3>& f,
        toolbox::Mesh<real_short,0>& fi,
        const std::array<int,D>& in,
        const std::array<int,D>& out
      );

  /// auxiliary functions to stagger e and b
  void stagger_x_eb(fields::YeeLattice& m);
  void stagger_y_eb(fields::YeeLattice& m);
  void stagger_z_eb(fields::YeeLattice& m);

  /// compute rho = div E
  void comp_rho(Tile<D>& tile);

  /// compute dE and dB
  void push_eb(Tile<D>& tile);

  /// compute and add jperp to E
  void add_jperp(Tile<D>& tile);

  /// update E and B
  void update_eb(Tile<D>& tile, real_short c1, real_short c2, real_short c3);

  /// set E_par = 0
  void remove_jpar(Tile<D>& tile);

  /// limit E < B
  void limit_e(Tile<D>& tile);

  /// copy Y^n to Y^n-1
  void copy_eb(Tile<D>& tile);

};



} // end of namespace

