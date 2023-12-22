#pragma once

#include "../../em-fields/tile.h"
#include "../../definitions.h"

namespace fields {

/// Rotating conductor
template<size_t D>
class Conductor
{

  public:

  Conductor() = default;

  // configuration parameters
  double radius = 10.0;
  double period = 0.0;
  double B0 = 1.0;  /// Initial magnetic field strength B_0
  double chi = 0.0; /// Obliquity angle between rotation axis and magnetic moment
  double phase = 0.0; /// rotator phase
  double cenx, ceny, cenz; // center of the sphere

  double delta = 1.0;

  // derived global quantities
  double angular_velocity;

  // additional switches for a possible polar cap 
  double radius_pc = 1.0; // polar cap radius
  double delta_pc  = 1.0; // polar smoothing


  // grid size (for nulling the sides to prevent periodic bc's)
  int Nx = 1;
  int Ny = 1;
  int Nz = 1;

  /// \brief interpolate electromagnetic fields to particle locations
  //void solve(fields::Tile<D>&  /*tile*/);

  double dipole(double x, double y, double z, int dim);

  void insert_em(fields::Tile<D>&  tile);

  void update_b(fields::Tile<D>&  tile);

  void update_e(fields::Tile<D>&  tile);

  void null_edges(fields::Tile<D>& tile, int mode);

};


} // end of namespace fields
