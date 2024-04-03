#pragma once

#include "core/ffe/tile.h"
#include "definitions.h"

namespace ffe {

/// Reduced FFE 4th order solver 
template<size_t D>
class FFE4
{
  public:

  int Nx;
  int Ny;
  int Nz;

  float eta = 1.0e-3; //resistivity for diffusion
  float reltime = 1.0; // e.b relaxation time (in units of dt)


  /// interpolated arrays
  toolbox::Mesh<float, 0> bxf;
  toolbox::Mesh<float, 0> byf;
  toolbox::Mesh<float, 0> bzf;

  toolbox::Mesh<float, 0> exf;
  toolbox::Mesh<float, 0> eyf;
  toolbox::Mesh<float, 0> ezf;

  toolbox::Mesh<float, 0> rhf;

  // extra arrays for jpar step
  toolbox::Mesh<float, 3> curlex;
  toolbox::Mesh<float, 3> curley;
  toolbox::Mesh<float, 3> curlez;

  toolbox::Mesh<float, 3> curlbx;
  toolbox::Mesh<float, 3> curlby;
  toolbox::Mesh<float, 3> curlbz;

  // extra arrays for jpar step
  toolbox::Mesh<float, 0> curlexf;
  toolbox::Mesh<float, 0> curleyf;
  toolbox::Mesh<float, 0> curlezf;
                             
  toolbox::Mesh<float, 0> curlbxf;
  toolbox::Mesh<float, 0> curlbyf;
  toolbox::Mesh<float, 0> curlbzf;


  FFE4(int Nx, int Ny, int Nz) :
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

  virtual ~FFE4() = default;

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

  /// remove jpar
  void remove_jpar(Tile<D>& tile);

  /// limit E < B
  void limit_e(Tile<D>& tile);

  ///  diffuse
  void add_diffusion(Tile<D>& tile);

};



} // end of namespace

