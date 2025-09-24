#pragma once

#include "core/emf/tile.h"
#include "definitions.h"
#include "tools/vector.h"


namespace emf {

using toolbox::Vec3;

/// Rotating conductor
template<size_t D>
class Conductor
{

  public:

  Conductor() = default;

  // configuration parameters
  float radius = 10.0; // stellar radius
  float period = 0.0;  // stellar rotation period
  float B0     = 1.0; // Initial magnetic field strength B_0
  float chi_mu = 0.0; // Obliquity angle of magnetic dipole from z-axis
  float chi_om = 0.0; // Obliquity angle of Omega vector from z-axis
  float phase_mu  = 0.0; // rotator phase of magnetic moment
  float phase_om  = 0.0; // rotator phase of rotation vector
  float cenx = 0, ceny = 0, cenz = 0; // center of the sphere

  float delta = 1.0;

  // derived global quantities
  float angular_velocity;

  // additional switches for a possible polar cap 
  float radius_pc = 1.0; // polar cap radius
  float delta_pc  = 1.0; // polar smoothing

  bool flat_surface = true; // use a flat (non-spherical/curved) surface for the star; makes BCs easier

  bool set_const_b = false; // oneD option to set background Bx to be constant (not a dipole)

  // grid size (for nulling the sides to prevent periodic bc's)
  int Nx = 1;
  int Ny = 1;
  int Nz = 1;

  /// \brief interpolate electromagnetic fields to particle locations
  //void solve(emf::Tile<D>&  /*tile*/);

  Vec3<float> dipole(Vec3<float>& xvec);

  void insert_em(emf::Tile<D>&  tile);

  void update_b(emf::Tile<D>&  tile);

  void update_e(emf::Tile<D>&  tile);

  void update_j(emf::Tile<D>&  tile);
};


} // end of namespace emf
