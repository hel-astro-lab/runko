#include <cmath> 
#include <cassert>
#include <algorithm>

#include "core/pic/boundaries/piston.h"
#include "tools/signum.h"

using std::min;
using std::max;

//TODO: FIXME direct copy from depositer/zigzag.c++
// use friending of the class or inherit this?
template<size_t D>
void pic::Piston<D>::zigzag(
    pic::Tile<D>& tile,
    /*old loc*/
    double x1glob, 
    double y1glob, 
    double z1glob, 
    /*new loc*/
    double x2glob, 
    double y2glob, 
    double z2glob, 
    double q)
{
  auto& gs = tile.get_grids();
  auto mins = tile.mins;

  // normalized location w.r.t. tile
  double x1 = D >= 1 ? x1glob - mins[0] : x1glob;
  double x2 = D >= 1 ? x2glob - mins[0] : x2glob;
  double y1 = D >= 2 ? y1glob - mins[1] : y1glob;
  double y2 = D >= 2 ? y2glob - mins[1] : y2glob;
  double z1 = D >= 3 ? z1glob - mins[2] : z1glob;
  double z2 = D >= 3 ? z2glob - mins[2] : z2glob;

  int i1  = D >= 1 ? floor( x1 ) : 0;
  int i2  = D >= 1 ? floor( x2 ) : 0;
  int j1  = D >= 2 ? floor( y1 ) : 0;
  int j2  = D >= 2 ? floor( y2 ) : 0;
  int k1  = D >= 3 ? floor( z1 ) : 0;
  int k2  = D >= 3 ? floor( z2 ) : 0;


  // relay point; +1 is equal to +\Delta x
  double xr = min( double(min(i1,i2)+1), max( double(max(i1,i2)), double(0.5*(x1+x2))) );
  double yr = min( double(min(j1,j2)+1), max( double(max(j1,j2)), double(0.5*(y1+y2))) );
  double zr = min( double(min(k1,k2)+1), max( double(max(k1,k2)), double(0.5*(z1+z2))) );


  // -q to include -j in the Ampere's equation
  double Fx1 = +q*(xr - x1);
  double Fy1 = +q*(yr - y1);
  double Fz1 = +q*(zr - z1);
  
  double Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
  double Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
  double Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

  double Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
  double Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  double Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;

  double Fx2 = +q*(x2-xr);
  double Fy2 = +q*(y2-yr);
  double Fz2 = +q*(z2-zr);

  // jx
  if (D >= 1) gs.jx(i1,  j1,   k1)   += Fx1 * (1.0-Wy1) * (1.0-Wz1);
  if (D >= 2) gs.jx(i1,  j1+1, k1)   += Fx1 * Wy1       * (1.0-Wz1);
  if (D >= 3) gs.jx(i1,  j1,   k1+1) += Fx1 * (1.0-Wy1) * Wz1;
  if (D >= 3) gs.jx(i1,  j1+1, k1+1) += Fx1 * Wy1     * Wz1;

  if (D >= 1) gs.jx(i2,  j2,   k2)   += Fx2 * (1.0-Wy2) * (1.0-Wz2);
  if (D >= 2) gs.jx(i2,  j2+1, k2)   += Fx2 * Wy2       * (1.0-Wz2);
  if (D >= 3) gs.jx(i2,  j2,   k2+1) += Fx2 * (1.0-Wy2) * Wz2;
  if (D >= 3) gs.jx(i2,  j2+1, k2+1) += Fx2 * Wy2     * Wz2;

  // jy
  if (D >= 1) gs.jy(i1,  j1,   k1)   += Fy1 * (1.0-Wx1) * (1.0-Wz1);
  if (D >= 2) gs.jy(i1+1,j1,   k1)   += Fy1 * Wx1       * (1.0-Wz1);
  if (D >= 3) gs.jy(i1  ,j1,   k1+1) += Fy1 * (1.0-Wx1) * Wz1;
  if (D >= 3) gs.jy(i1+1,j1,   k1+1) += Fy1 * Wx1     * Wz1;
  
  if (D >= 1) gs.jy(i2,  j2,   k2)   += Fy2 * (1.0-Wx2) * (1.0-Wz2);
  if (D >= 2) gs.jy(i2+1,j2,   k2)   += Fy2 * Wx2       * (1.0-Wz2);
  if (D >= 3) gs.jy(i2,  j2,   k2+1) += Fy2 * (1.0-Wx2) * Wz2;
  if (D >= 3) gs.jy(i2+1,j2,   k2+1) += Fy2 * Wx2     * Wz2;
                        
  // jz
  if (D >= 1) gs.jz(i1,  j1,   k1)   += Fz1 * (1.0-Wx1) * (1.0-Wy1);
  if (D >= 1) gs.jz(i1+1,j1,   k1)   += Fz1 * Wx1       * (1.0-Wy1);
  if (D >= 2) gs.jz(i1,  j1+1, k1)   += Fz1 * (1.0-Wx1) * Wy1;
  if (D >= 2) gs.jz(i1+1,j1+1, k1)   += Fz1 * Wx1       * Wy1;

  if (D >= 2) gs.jz(i2,  j2,   k2)   += Fz2 * (1.0-Wx2) * (1.0-Wy2);
  if (D >= 2) gs.jz(i2+1,j2,   k2)   += Fz2 * Wx2       * (1.0-Wy2);
  if (D >= 2) gs.jz(i2,  j2+1, k2)   += Fz2 * (1.0-Wx2) * Wy2;
  if (D >= 2) gs.jz(i2+1,j2+1, k2)   += Fz2 * Wx2       * Wy2;

}

template<size_t D>
void pic::Piston<D>::solve(
    pic::Tile<D>& tile)
{

  // outflowing particles
  std::vector<int> to_be_deleted;

  // skip if piston head is not inside tile boundaries
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  if(!(mins[0] <= walloc && walloc <= maxs[0])) return;


  for(auto&& container : tile.containers) {
    to_be_deleted.clear();

    const double c = tile.cfl;
    const double q = container.q;

    for(size_t n=0; n<container.size(); n++) {

      // left side of the wall boundary
      if( container.loc(0,n) < walloc) {

        //prtcl step and velocity at updated time
        double x1 = container.loc(0,n);
        double y1 = container.loc(1,n);
        double z1 = container.loc(2,n);
                                                     
        double u1 = container.vel(0,n);
        double v1 = container.vel(1,n);
        double w1 = container.vel(2,n);

        // unwind wall location one time step backwards
        double walloc0 = walloc - betawall*c;

        // unwind prtcl location one time step backwards
        double gamma = sqrt(1.0 + u1*u1 + v1*v1 + w1*w1);
        // NOTE: this could be approximated by only treating the refl. component
        double x0 = x1 - c*u1/gamma;
        double y0 = y1 - c*v1/gamma;
        double z0 = z1 - c*w1/gamma;

        // check if this is reflected particle; else remove 
        // equals to right boundary outflow as they wrap particles to here
        if(walloc0 - x0 > c) {
          to_be_deleted.push_back(n);
          continue;
        }

        //--------------------------------------------------
        // compute crossing point
        
        // time step fraction spent traveling towards the reflector 
        double dt = std::abs((x0-walloc0)/(betawall*c - c*u1/gamma + EPS));

        // skip this particle since it is further away from the wall than one time step
        if(dt > 1.0) {
          to_be_deleted.push_back(n);
          continue; 
        }

        // NOTE: this could be approximated by only treating the refl. component
        double xcolis = x0 + c*dt*u1/gamma;
        double ycolis = y0 + c*dt*v1/gamma;
        double zcolis = z0 + c*dt*w1/gamma;
          
        // deposit current from beginning upto intersection w/ wall
        zigzag(tile, x0, y0, z0, xcolis, ycolis, zcolis, q);

        //--------------------------------------------------
        // perform reflection

        // rewrite particle momentum as it gets a kick from the wall
        //u1 = -u1; // simple reflection w/o wall movement
        u1 = gammawall*gammawall*gamma*(2.*betawall - u1/gamma*(1.0+betawall*betawall));

        gamma = sqrt(1.0 + u1*u1 + v1*v1 + w1*w1);

        // NOTE: we use x1-x0 instead of c*vx/gam to get the old vel before kick
        // this is recomputed because prtcl gamma might change if the wall is moving
        dt = std::min(1.0, std::abs((x1 - xcolis)/(x1-x0 + EPS)) );

        // move particle from the location of the wall with new velocity
        double xnew = xcolis + c*dt*u1/gamma;
        double ynew = ycolis + c*dt*v1/gamma;
        double znew = zcolis + c*dt*w1/gamma;


        // clean up the part of trajectory behind the wall that will be added by the 
        // deposition routine that unwinds the particle location by a full time step.
        //
        // - loc-c*u/gamma is the unwinded starting location, x_n performed by the normal zigzag
        // - xcolis is the collision point.
        // - injecting opposite charge to this region to cancel out normal zigzag current
        // - rest of the (from wall to new pos) is done by the normal deposit call
        zigzag(tile, 
            xnew - c*u1/gamma,
            ynew - c*v1/gamma,
            znew - c*w1/gamma,
            xcolis, ycolis, zcolis,
            -q);

        // lastly; store particle back to the container
        container.loc(0,n) = static_cast<float>( xnew );
        container.loc(1,n) = static_cast<float>( ynew );
        container.loc(2,n) = static_cast<float>( znew );

        container.vel(0,n) = static_cast<float>( u1 );
        //container.vel(1,n) = static_cast<float>( v1 );
        //container.vel(2,n) = static_cast<float>( w1 );
      }
    }

    // process outflown particles
    container.delete_particles(to_be_deleted);

  } // end of loop over species
}


template<>
void pic::Piston<1>::field_bc(
    pic::Tile<1>& tile)
{
  int k = 0; // collapse third D
  int j = 0; // collapse second D

  // skip if piston head is not inside tile boundaries
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  // make left side of piston conductor
  if(walloc < maxs[0]) {
    auto& gs = tile.get_grids();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0]);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero
    for(int i=-3; i<=iw; i++) {

        // transverse components of electric field to zero (only parallel comp allowed)
        gs.ey(i,j,k) = 0.0;
        gs.ez(i,j,k) = 0.0;

        // clean all current behind piston head
        // NOTE: not needed since we null the current inside reflection procedure
        //gs.jx(i,j,k) = 0.0;
        //gs.jy(i,j,k) = 0.0;
        //gs.jz(i,j,k) = 0.0;
    }
  }

}


template<>
void pic::Piston<2>::field_bc(
    pic::Tile<2>& tile)
{
  int k = 0; // collapse third D

  // skip if piston head is not inside tile boundaries
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  // make left side of piston conductor
  if(walloc < maxs[0]) {
    auto& gs = tile.get_grids();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0]);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) {
      for(int i=-3; i<=iw; i++) {

        // transverse components of electric field to zero (only parallel comp allowed)
        gs.ey(i,j,k) = 0.0;
        gs.ez(i,j,k) = 0.0;

        // clean all current behind piston head
        // NOTE: not needed since we null the current inside reflection procedure
        //gs.jx(i,j,k) = 0.0;
        //gs.jy(i,j,k) = 0.0;
        //gs.jz(i,j,k) = 0.0;
      }
    }
  }

}


template<>
void pic::Piston<3>::field_bc(
    pic::Tile<3>& tile)
{
  // skip if piston head is not inside tile boundaries
  auto mins = tile.mins;
  auto maxs = tile.maxs;

  // make left side of piston conductor
  if(walloc < maxs[0]) {
    auto& gs = tile.get_grids();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0]);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero to make this conductor
    for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
    for(int i=-3; i<=iw; i++) {

      // transverse components of electric field to zero (only parallel comp allowed)
      gs.ey(i,j,k) = 0.0;
      gs.ez(i,j,k) = 0.0;

      // clean all current behind piston head
      //gs.jx(i,j,k) = 0.0;
      //gs.jy(i,j,k) = 0.0;
      //gs.jz(i,j,k) = 0.0;
    }

  } // end if inside piston

}


//--------------------------------------------------
// explicit template instantiation

template class pic::Piston<1>; // 1D3V
template class pic::Piston<2>; // 2D3V
template class pic::Piston<3>; // 3D3V

