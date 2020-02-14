#pragma once

#include "../tile.h"
#include "../../definitions.h"

namespace ffe {

/// Combined RK3 reduced FFE 2nd order solver 
template<size_t D>
class rFFE2
{

  int Nx;
  int Ny;
  int Nz;

  /// interpolated arrays
  toolbox::Mesh<real_short, 0> ibx;
  toolbox::Mesh<real_short, 0> iby;
  toolbox::Mesh<real_short, 0> ibz;

  toolbox::Mesh<real_short, 0> iex;
  toolbox::Mesh<real_short, 0> iey;
  toolbox::Mesh<real_short, 0> iez;

  toolbox::Mesh<real_short, 0> irh;


  rFFE2(int Nx, int Ny, int Nz) :
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz)
    exf(Nx, Ny, Nz),
    eyf(Nx, Ny, Nz),
    ezf(Nx, Ny, Nz),
    rhf(Nx, Ny, Nz)
  {};

  virtual ~rFFE2() = default;


  // interpolation routine
  void interpolate( 
        toolbox::Mesh<real_short,3>& f,
        std::array<int,D>& in,
        std::array<int,D>& out
      );

  /// compute dE and dB
  void push_eb(Tile<D>& tile);

  /// compute and add jperp to E
  void add_jperp(Tile<D>& tile);

  /// update E and B
  void update_eb(Tile<D>& tile);

  /// set E_par = 0
  void remove_epar(Tile<D>& tile);

  /// limit E < B
  void limit_e(Tile<D>& tile);

};


} // end of namespace

