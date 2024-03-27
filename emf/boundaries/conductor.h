#pragma once

#include "emf/tile.h"
#include "definitions.h"
#include "tools/vector.h"


namespace emf {

using toolbox::Vec3;


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
  float sx,sy,sz;
  float cx, cy, cz, r;

  public:

  StaggeredSphericalField( 
      float sx, float sy, float sz,
      float cenx, float ceny, float cenz, float radius) 
    : sx(sx), sy(sy), sz(sz),
      cx(cenx), cy(ceny), cz(cenz), r(radius)
  { }

  // NOTE there is a flip of staggering direction for negative cartesian coordinates
  //      not really sure why, but it is needed to get a balanced configuration
  //inline float x(float i) { return i > cx ? (i + 0.5*sx - cx)/r : (i - 0.5*sx - cx)/r; }
  //inline float y(float j) { return j > cy ? (j + 0.5*sy - cy)/r : (j - 0.5*sy - cy)/r; }
  //inline float z(float k) { return k > cz ? (k + 0.5*sz - cz)/r : (k - 0.5*sz - cz)/r; }

  inline float x(float i) { return (i + 0.5*sx - cx)/r; }
  inline float y(float j) { return (j + 0.5*sy - cy)/r; }
  inline float z(float k) { return (k + 0.5*sz - cz)/r; }

  // TODO template dependency complicates the whole code too much
  //      feeding D in as a parameter for simplicity
  inline toolbox::Vec3<float> vec(float i, float j, float k, size_t D) 
  { 
    if(D == 1) return Vec3<float>( x(i),  0.0,  0.0 );
    if(D == 2) return Vec3<float>( x(i), y(j),  0.0 );
    if(D == 3) return Vec3<float>( x(i), y(j), z(k) );
    assert(false);
  }
};


class StaggeredSphericalCoordinates
{
  float cx, cy, cz, r;

  public:

  StaggeredSphericalCoordinates( float cenx, float ceny, float cenz, float radius) 
      : cx(cenx), cy(ceny), cz(cenz), r(radius)
  {}

  inline StaggeredSphericalField mid(){ return StaggeredSphericalField(1., 1., 1., cx,cy,cz,r); }

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

/*
   smooth ramp; half-way is at r = r0; 
   r0 +- 2*delta is roughly the asymptotic regime of the function as shown below
  
   1 | ----
     |      \
     |       \
   0 |         ------
     -------|----------->
            r0
*/
inline float_m shape(float_m r, float_m r0, float_m delta) 
{ 
  return 0.5 * (1.0 - tanh( (r - r0) / delta )); 
}


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

  bool flat_surface = true;

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

};


} // end of namespace emf
