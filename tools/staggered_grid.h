#pragma once

#include "definitions.h"
#include "tools/vector.h"

using toolbox::Vec3;

namespace toolbox {

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
inline float shape(float r, float r0, float delta) 
{ 
  return 0.5 * (1.0 - tanh( (r - r0) / delta )); 
}


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



} // namespace
