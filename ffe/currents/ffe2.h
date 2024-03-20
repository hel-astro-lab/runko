#pragma once

#include "../tile.h"
#include "../../definitions.h"

namespace ffe {

/// Reduced FFE 4th order solver 
template<size_t D>
class FFE2
{
  public:

  int Nx;
  int Ny;
  int Nz;

  float_m eta = 1.0e-3; //resistivity for diffusion
  float_m reltime = 1.0; // e.b relaxation time (in units of dt)

  /// interpolated arrays
  toolbox::Mesh<float_m, 0> bxf;
  toolbox::Mesh<float_m, 0> byf;
  toolbox::Mesh<float_m, 0> bzf;

  toolbox::Mesh<float_m, 0> exf;
  toolbox::Mesh<float_m, 0> eyf;
  toolbox::Mesh<float_m, 0> ezf;

  toolbox::Mesh<float_m, 0> rhf;

  // extra arrays for jpar step
  toolbox::Mesh<float_m, 3> curlex;
  toolbox::Mesh<float_m, 3> curley;
  toolbox::Mesh<float_m, 3> curlez;

  toolbox::Mesh<float_m, 3> curlbx;
  toolbox::Mesh<float_m, 3> curlby;
  toolbox::Mesh<float_m, 3> curlbz;

  // extra arrays for jpar step
  toolbox::Mesh<float_m, 0> curlexf;
  toolbox::Mesh<float_m, 0> curleyf;
  toolbox::Mesh<float_m, 0> curlezf;
                             
  toolbox::Mesh<float_m, 0> curlbxf;
  toolbox::Mesh<float_m, 0> curlbyf;
  toolbox::Mesh<float_m, 0> curlbzf;


  FFE2(int Nx, int Ny, int Nz) :
    Nx(Nx), Ny(Ny), Nz(Nz),
    bxf(Nx, Ny, Nz),
    byf(Nx, Ny, Nz),
    bzf(Nx, Ny, Nz),
    exf(Nx, Ny, Nz),
    eyf(Nx, Ny, Nz),
    ezf(Nx, Ny, Nz),
    rhf(Nx, Ny, Nz),
    curlex(Nx, Ny, Nz),
    curley(Nx, Ny, Nz),
    curlez(Nx, Ny, Nz),
    curlbx(Nx, Ny, Nz),
    curlby(Nx, Ny, Nz),
    curlbz(Nx, Ny, Nz),
    curlexf(Nx, Ny, Nz),
    curleyf(Nx, Ny, Nz),
    curlezf(Nx, Ny, Nz),
    curlbxf(Nx, Ny, Nz),
    curlbyf(Nx, Ny, Nz),
    curlbzf(Nx, Ny, Nz)
  {};

  virtual ~FFE2() = default;

  // interpolation routine
  void interpolate( 
        toolbox::Mesh<float_m,3>& f,
        toolbox::Mesh<float_m,0>& fi,
        const std::array<int,D>& in,
        const std::array<int,D>& out
      );

  /// auxiliary functions to stagger e and b
  void stagger_x_eb(emf::YeeLattice& m);
  void stagger_y_eb(emf::YeeLattice& m);
  void stagger_z_eb(emf::YeeLattice& m);

  void stagger_x_curl();
  void stagger_y_curl();
  void stagger_z_curl();

  /// compute rho = div E
  void comp_rho(Tile<D>& tile);

  /// compute dE and dB
  void push_eb(Tile<D>& tile);

  /// compute and add jperp to E
  void add_jperp(Tile<D>& tile);

  /// compute jpar and add to E
  void add_jpar(Tile<D>& tile);

  /// limit E < B
  void limit_e(Tile<D>& tile);

  /// diffusion
  void add_diffusion(Tile<D>& tile);

};



} // end of namespace

