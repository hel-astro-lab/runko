#include "conductor.h"

#include <cmath> 
#include <cassert>

using std::min;
using std::max;


// General dipole formula in cartesian coordinates
template<>
double fields::Conductor<3>::dipole(
        double x,
        double y,
        double z,
        int dim) {

  //non-rotating magnetic moment vector components
  double p1 = sin(chi), p2 = 0.0, p3 = cos(chi);

  // final rotating magnetic moment with phase included
  double mux = p1*cos(phase) - p2*sin(phase);
  double muy = p1*sin(phase) + p2*cos(phase);
  double muz = p3;

  double rad = std::sqrt(x*x + y*y + z*z);

  // mu . r
  double mudotr = mux*x + muy*y + muz*z;

   if     (dim == 0) return 3.0*x*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS); //x
   else if(dim == 1) return 3.0*y*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS); //y
   else if(dim == 2) return 3.0*z*mudotr/( pow(rad,5) + EPS) - muz/(pow(rad,3) + EPS); //z
   return 0.0;
}


//real_short S(real_short r, real_short r0, real_short delta) { return 0.5 * (1.0 - tanh((r - r0) / delta)); }



template<>
void fields::Conductor<3>::insert_em(
    fields::Tile<3>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  real_short bxd, byd, bzd, exd, eyd, ezd;
  real_short iglob, jglob, kglob;
  real_short xr,yr,zr;

  // set transverse directions to zero to make this conductor
  for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
  for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
  for(int i=-3; i<static_cast<int>(tile.mesh_lengths[0])+3; i++) {

    //-----------
    // coordinates
    // TODO staggering
      
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];

    // x coord: 1,0,0
    // y coord: 0,1,0
    // z coord: 0,0,1
    
    xr = (iglob - cenx)/radius;
    yr = (jglob - ceny)/radius;
    zr = (kglob - cenz)/radius;

    //-----------
    // magnetic field
    
    bxd = B0*dipole(xr,yr,zr,0);
    byd = B0*dipole(xr,yr,zr,1);
    bzd = B0*dipole(xr,yr,zr,2);

    yee.bx(i,j,k) = bxd;
    yee.by(i,j,k) = byd;
    yee.bz(i,j,k) = bzd;

    //-----------
    // electric field

    exd = 0.0;
    eyd = 0.0;
    ezd = 0.0;

    yee.ex(i,j,k) = exd;
    yee.ey(i,j,k) = eyd;
    yee.ez(i,j,k) = ezd;

    //yee.jx(i,j,k) = 0.0;
    //yee.jy(i,j,k) = 0.0;
    //yee.jz(i,j,k) = 0.0;
  }

}

template<>
void fields::Conductor<3>::update_b(
    fields::Tile<3>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  // set transverse directions to zero to make this conductor
  for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
  for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
  for(int i=-3; i<static_cast<int>(tile.mesh_lengths[0])+3; i++) {

    yee.bx(i,j,k) = 0.0;
    yee.by(i,j,k) = 0.0;
    yee.bz(i,j,k) = 0.0;

    // clean all current behind piston head
    //yee.jx(i,j,k) = 0.0;
    //yee.jy(i,j,k) = 0.0;
    //yee.jz(i,j,k) = 0.0;
  }

}


template<>
void fields::Conductor<3>::update_e(
    fields::Tile<3>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  // set transverse directions to zero to make this conductor
  for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
  for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
  for(int i=-3; i<static_cast<int>(tile.mesh_lengths[0])+3; i++) {
    yee.ex(i,j,k) = 0.0;
    yee.ey(i,j,k) = 0.0;
    yee.ez(i,j,k) = 0.0;

    // clean all current behind piston head
    yee.jx(i,j,k) = 0.0;
    yee.jy(i,j,k) = 0.0;
    yee.jz(i,j,k) = 0.0;
  }

}


//--------------------------------------------------
// explicit template instantiation

//template class fields::Conductor<2>; // 2D
template class fields::Conductor<3>; // 3D

