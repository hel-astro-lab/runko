#include "piston.h"
#include "../../tools/signum.h"

#include <cmath> 
#include <cassert>
#include <algorithm>

using std::min;
using std::max;

//TODO: FIXME direct copy from depositer/zigzag.c++
// use friending of the class or inherit this?
template<size_t D>
void pic::Piston<D>::zigzag(
    pic::Tile<D>& tile,
    /*old loc*/
    real_long x1glob, 
    real_long y1glob, 
    real_long z1glob, 
    /*new loc*/
    real_long x2glob, 
    real_long y2glob, 
    real_long z2glob, 
    real_long q)
{
  auto& yee = tile.get_yee();
  auto mins = tile.mins;

  // normalized location w.r.t. tile
  real_long x1 = D >= 1 ? x1glob - mins[0] : x1glob;
  real_long x2 = D >= 1 ? x2glob - mins[0] : x2glob;
  real_long y1 = D >= 2 ? y1glob - mins[1] : y1glob;
  real_long y2 = D >= 2 ? y2glob - mins[1] : y2glob;
  real_long z1 = D >= 3 ? z1glob - mins[2] : z1glob;
  real_long z2 = D >= 3 ? z2glob - mins[2] : z2glob;

  int i1  = D >= 1 ? static_cast<int>( floor(x1) ) : 0;
  int i2  = D >= 1 ? static_cast<int>( floor(x2) ) : 0;
  int j1  = D >= 2 ? static_cast<int>( floor(y1) ) : 0;
  int j2  = D >= 2 ? static_cast<int>( floor(y2) ) : 0;
  int k1  = D >= 3 ? static_cast<int>( floor(z1) ) : 0;
  int k2  = D >= 3 ? static_cast<int>( floor(z2) ) : 0;


  // relay point; +1 is equal to +\Delta x
  real_long xr = min( real_long(min(i1,i2)+1), max( real_long(max(i1,i2)), real_long(0.5*(x1+x2)) ) );
  real_long yr = min( real_long(min(j1,j2)+1), max( real_long(max(j1,j2)), real_long(0.5*(y1+y2)) ) );
  real_long zr = min( real_long(min(k1,k2)+1), max( real_long(max(k1,k2)), real_long(0.5*(z1+z2)) ) );

  // -q to include -j in the Ampere's equation
  real_long Fx1 = +q*(xr - x1);
  real_long Fy1 = +q*(yr - y1);
  real_long Fz1 = +q*(zr - z1);
  
  real_long Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
  real_long Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
  real_long Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

  real_long Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
  real_long Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  real_long Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;

  // TODO: changed sign of q
  real_long Fx2 = +q*(x2-xr);
  real_long Fy2 = +q*(y2-yr);
  real_long Fz2 = +q*(z2-zr);

  // jx
  if (D >= 1) yee.jx(i1,  j1,   k1)   += Fx1 * (1.0-Wy1) * (1.0-Wz1);
  if (D >= 2) yee.jx(i1,  j1+1, k1)   += Fx1 * Wy1       * (1.0-Wz1);
  if (D >= 3) yee.jx(i1,  j1,   k1+1) += Fx1 * (1.0-Wy1) * Wz1;
  if (D >= 3) yee.jx(i1,  j1+1, k1+1) += Fx1 * Wy1       * Wz1;

  if (D >= 1) yee.jx(i2,  j2,   k2)   += Fx2 * (1.0-Wy2) * (1.0-Wz2);
  if (D >= 2) yee.jx(i2,  j2+1, k2)   += Fx2 * Wy2       * (1.0-Wz2);
  if (D >= 3) yee.jx(i2,  j2,   k2+1) += Fx2 * (1.0-Wy2) * Wz2;
  if (D >= 3) yee.jx(i2,  j2+1, k2+1) += Fx2 * Wy2       * Wz2;

  // jy
  if (D >= 1) yee.jy(i1,  j1,   k1)   += Fy1 * (1.0-Wx1) * (1.0-Wz1);
  if (D >= 2) yee.jy(i1+1,j1,   k1)   += Fy1 * Wx1       * (1.0-Wz1);
  if (D >= 3) yee.jy(i1  ,j1,   k1+1) += Fy1 * (1.0-Wx1) * Wz1;
  if (D >= 3) yee.jy(i1+1,j1,   k1+1) += Fy1 * Wx1       * Wz1;
  
  if (D >= 1) yee.jy(i2,  j2,   k2)   += Fy2 * (1.0-Wx2) * (1.0-Wz2);
  if (D >= 2) yee.jy(i2+1,j2,   k2)   += Fy2 * Wx2       * (1.0-Wz2);
  if (D >= 3) yee.jy(i2,  j2,   k2+1) += Fy2 * (1.0-Wx2) * Wz2;
  if (D >= 3) yee.jy(i2+1,j2,   k2+1) += Fy2 * Wx2       * Wz2;
                        
  // jz
  if (D >= 1) yee.jz(i1,  j1,   k1)   += Fz1 * (1.0-Wx1) * (1.0-Wy1);
  if (D >= 1) yee.jz(i1+1,j1,   k1)   += Fz1 * Wx1       * (1.0-Wy1);
  if (D >= 2) yee.jz(i1,  j1+1, k1)   += Fz1 * (1.0-Wx1) * Wy1;
  if (D >= 2) yee.jz(i1+1,j1+1, k1)   += Fz1 * Wx1       * Wy1;

  if (D >= 1) yee.jz(i2,  j2,   k2)   += Fz2 * (1.0-Wx2) * (1.0-Wy2);
  if (D >= 1) yee.jz(i2+1,j2,   k2)   += Fz2 * Wx2       * (1.0-Wy2);
  if (D >= 2) yee.jz(i2,  j2+1, k2)   += Fz2 * (1.0-Wx2) * Wy2;
  if (D >= 2) yee.jz(i2+1,j2+1, k2)   += Fz2 * Wx2       * Wy2;

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

  real_long wshf = 0.0; // staggering shift for the prtcl wall
  if(!(mins[0] <= walloc+wshf && walloc+wshf <= maxs[0])) return;


  for(auto&& container : tile.containers) {
    to_be_deleted.clear();

    const real_long c = tile.cfl;
    const real_long q = container.q;

    for(int n=0; n<static_cast<int>(container.size()); n++) {

      // left side of the wall boundary
      if( container.loc(0,n) < walloc+wshf) {

        //prtcl step and velocity at updated time
        real_long x1 = container.loc(0,n);
        real_long y1 = container.loc(1,n);
        real_long z1 = container.loc(2,n);
                                                     
        real_long u1 = container.vel(0,n);
        real_long v1 = container.vel(1,n);
        real_long w1 = container.vel(2,n);

        // unwind wall location one time step backwards
        real_long walloc0 = walloc +wshf - betawall*c;

        // unwind prtcl location one time step backwards
        real_long gamma = sqrt(1.0 + u1*u1 + v1*v1 + w1*w1);
        // NOTE: this could be approximated by only treating the refl. component
        real_long x0 = x1 - c*u1/gamma;
        real_long y0 = y1 - c*v1/gamma;
        real_long z0 = z1 - c*w1/gamma;

        // check if this is reflected particle; else remove 
        // equals to right boundary outflow as it wraps particles here
        if(walloc0 - x0 > c) {
          to_be_deleted.push_back(n);
          continue;
        }

        //--------------------------------------------------
        // compute crossing point
        
        // time step fraction spent traveling towards the reflector 
        //real_long dt = std::abs((x0-walloc0)/(betawall*c - c*u1/gamma));
        real_long wvel = std::max(EPS, std::abs(betawall*c - c*u1/gamma)); 
        real_long dt = std::abs(x0-walloc0)/wvel;
        dt = std::min(1.0, std::abs(x0-walloc0)/wvel);

        // NOTE: this could be approximated by only treating the refl. component
        real_long xcolis = x0 + c*dt*u1/gamma;
        real_long ycolis = y0 + c*dt*v1/gamma;
        real_long zcolis = z0 + c*dt*w1/gamma;
          
        // deposit current from beginning upto intersection w/ wall
        zigzag(tile, x0, y0, z0, xcolis, ycolis, zcolis, +q);

        //--------------------------------------------------
        // perform reflection

        // rewrite particle momentum as it gets a kick from the wall
        //u1 = -u1; // simple reflection w/o wall movement
        u1 = gammawall*gammawall*gamma*(2.*betawall - u1/gamma*(1.0+betawall*betawall));
        gamma = sqrt(1.0 + u1*u1 + v1*v1 + w1*w1);

        // NOTE: we use x1-x0 instead of c*vx/gam to get the old vel before kick
        // this is recomputed because prtcl gamma might change if the wall is moving
        //dt = std::min(1.0, std::max(EPS, std::abs((x1 - xcolis)/(x1-x0))) );
        wvel = std::max(EPS, std::abs(x1-x0)); 
        dt = std::abs( (x1-xcolis)/wvel );
        dt = std::min(1.0, dt);

        // move particle from the location of the wall with new velocity
        real_long xnew = xcolis + c*dt*u1/gamma;
        real_long ynew = ycolis + c*dt*v1/gamma;
        real_long znew = zcolis + c*dt*w1/gamma;

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
        container.loc(0,n) = xnew;
        container.loc(1,n) = ynew;
        container.loc(2,n) = znew;

        container.vel(0,n) = u1;
        //container.vel(1,n) = v1;
        //container.vel(2,n) = w1;
      }
    }

    // process outflown particles
    container.delete_particles(to_be_deleted);

  } // end of loop over species
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
  if(walloc -10.0 < maxs[0]) {
    auto& yee = tile.get_yee();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0] - 10.0);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) {
      for(int i=-3; i<=iw; i++) {

        // transverse components of electric field to zero (only parallel comp allowed)
        yee.ey(i,j,k) = 0.0;
        yee.ez(i,j,k) = 0.0;

        // clean all current behind piston head
        //yee.ex(i,j,k) = 0.0;
        //yee.jx(i,j,k) = 0.0;
        //yee.jy(i,j,k) = 0.0;
        //yee.jz(i,j,k) = 0.0;
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
  //auto maxs = tile.maxs;
      
  // NOTE: fields are set -10 cells behind prtcl reflection boundary
  double walloc0 = walloc - 10.0; // location of conductor 

  // make left side of piston conductor
  if(mins[0] < walloc0 -3) {
    auto& yee = tile.get_yee();

    // wall location 
    auto iw = static_cast<int>( floor(walloc0 - mins[0]) ); // to tile units
    if(iw > static_cast<int>(tile.mesh_lengths[0]+3)) iw = tile.mesh_lengths[0]+3;


    // set transverse directions to zero to make this conductor
    for(int i=-3; i<=iw; i++) {

    //std::cout << "piston field bc: ix" << i + tile.mins[0] << " iw" << iw << "vals" << eywall << " " << ezwall << "\n";

    for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) {
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) {

      // transverse components of electric field to zero (only parallel comp allowed)
      //yee.ey(i,j,k) = 0.0; // JN: testing what happens 4th of mar
      //yee.ez(i,j,k) = 0.0;
        
      // JN: TODO test 2
      //yee.bx(i,j,k) = bxwall;
      //yee.by(i,j,k) = bywall; 
      //yee.bz(i,j,k) = bzwall;
      //
      //yee.ex(i,j,k) = exwall;
      //yee.ey(i,j,k) = eywall; 
      //yee.ez(i,j,k) = ezwall;

      // JN: TODO test 3; DONE: nulling bx makes it even worse
      //yee.bx(i,j,k) = bxwall; 
      //yee.ey(i,j,k) = eywall; 
      //yee.ez(i,j,k) = ezwall;

      // JN: TODO test 4
      //yee.ey(i,  j,k) = eywall; 
      //yee.ez(i+1,j,k) = ezwall;

      
      // JN: TODO test 5
      yee.ey(i,j,k) = eywall; 
      yee.ez(i,j,k) = ezwall;

      // JN: TODO test 5
      // prevent time evolution of B field
      //yee.bx(i,j,k) = bxwall;
      //yee.by(i+1,j,k) = bywall; 
      //yee.bz(i+1,j,k) = bzwall;
      //
      // kill surface electric fields on conductor
      //yee.ey(i,j,k) = eywall; 
      //yee.ez(i,j,k) = ezwall;

      // clean all current behind piston head
      //yee.jx(i,j,k) = 0.0;
      //yee.jy(i,j,k) = 0.0;
      //yee.jz(i,j,k) = 0.0;
    }}}

  } // end if inside piston

}


//--------------------------------------------------
// explicit template instantiation

//template class pic::Piston<1>; // 1D3V
template class pic::Piston<2>; // 2D3V
template class pic::Piston<3>; // 3D3V

