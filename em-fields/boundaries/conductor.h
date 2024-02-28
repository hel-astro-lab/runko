#pragma once

#include "../../em-fields/tile.h"
#include "../../definitions.h"
#include "../../tools/vector.h"

// Helper class to automate the management of a staggered Yee lattice
// 
// Returns a field object that has the correct staggering. 
//
// Call syntax is:
//
//    StaggeredCoordinates coord();
//    auto x = coord.ex().x(i);

class StaggeredSphericalField
{
  double sx,sy,sz;
  double cx, cy, cz, r;

  public:

  StaggeredSphericalField( 
      double sx, double sy, double sz,
      double cenx, double ceny, double cenz, double radius) 
    : sx(sx), sy(sy), sz(sz),
      cx(cenx), cy(ceny), cz(cenz), r(radius)
  { }

  // NOTE there is a flip of staggering direction for negative cartesian coordinates
  //      not really sure why, but it is needed to get a balanced configuration
  inline double x(double i) { return i > cx ? (i + sx - cx)/r : (i - sx - cx)/r; }
  inline double y(double j) { return j > cy ? (j + sy - cy)/r : (j - sy - cy)/r; }
  inline double z(double k) { return k > cz ? (k + sz - cz)/r : (k - sz - cz)/r; }

  // TODO template dependency complicates the whole code too much
  //inline toolbox::Vec3<double> vec(double i, double j, double k) 
  //{ 
  //  if(D == 1) return ret( x(i),  0.0,  0.0 );
  //  if(D == 2) return ret( x(i), y(j),  0.0 );
  //  if(D == 3) return ret( x(i), y(j), z(k) );
  //}

  //inline double x(double i) { return (i + sx - cx)/r; }
  //inline double y(double j) { return (j + sy - cy)/r; }
  //inline double z(double k) { return (k + sz - cz)/r; }

  //void _dummy(std::array<int, D> in) = 0;
};


class StaggeredSphericalCoordinates
{
  double cx, cy, cz, r;

  public:

  StaggeredSphericalCoordinates( double cenx, double ceny, double cenz, double radius) 
      : cx(cenx), cy(ceny), cz(cenz), r(radius)
  {}

  inline StaggeredSphericalField rh() { return StaggeredSphericalField(0., 0., 0., cx,cy,cz,r); }

  inline StaggeredSphericalField ex() { return StaggeredSphericalField(1., 0., 0., cx,cy,cz,r); }
  inline StaggeredSphericalField ey() { return StaggeredSphericalField(0., 1., 0., cx,cy,cz,r); }
  inline StaggeredSphericalField ez() { return StaggeredSphericalField(1., 0., 1., cx,cy,cz,r); }

  inline StaggeredSphericalField jx() { return StaggeredSphericalField(1., 0., 0., cx,cy,cz,r); }
  inline StaggeredSphericalField jy() { return StaggeredSphericalField(0., 1., 0., cx,cy,cz,r); }
  inline StaggeredSphericalField jz() { return StaggeredSphericalField(1., 0., 1., cx,cy,cz,r); }

  inline StaggeredSphericalField bx() { return StaggeredSphericalField(0., 1., 1., cx,cy,cz,r); }
  inline StaggeredSphericalField by() { return StaggeredSphericalField(1., 0., 1., cx,cy,cz,r); }
  inline StaggeredSphericalField bz() { return StaggeredSphericalField(1., 1., 0., cx,cy,cz,r); }

};


// smooth ramp; half-way is at r = r0; 
// r0 +- 2*delta is roughly the asymptotic regime of the function as shown below
//
// 1 | ----
//   |      \
//   |       \
// 0 |         ------
//   -------|----------->
//          r0
inline float_m shape(float_m r, float_m r0, float_m delta) 
{ 
  return 0.5 * (1.0 - tanh( (r - r0) / delta )); 
}


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
  double B0     = 1.0;  /// Initial magnetic field strength B_0
  double chi_mu = 0.0; /// Obliquity angle of magnetic dipole from z-axis
  double chi_om = 0.0; /// Obliquity angle of Omega vector from z-axis
  double phase  = 0.0; /// rotator phase
  double cenx = 0, ceny = 0, cenz = 0; // center of the sphere

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
