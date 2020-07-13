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
    real_long x2glob, 
    real_long y2glob, 
    real_long z2glob, 
    real_long x1glob, 
    real_long y1glob, 
    real_long z1glob, 
    real_long q)
{
  int i1,i2,j1,j2,k1,k2;

  real_long x1,x2,y1,y2,z1,z2;
  real_long xr,yr,zr;

  real_long Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
  real_long Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;

  auto& yee = tile.get_yee();
  auto mins = tile.mins;

  // normalized location w.r.t. tile
  x1 = D >= 1 ? x1glob - mins[0] : x1glob;
  x2 = D >= 1 ? x2glob - mins[0] : x2glob;
  y1 = D >= 2 ? y1glob - mins[1] : y1glob;
  y2 = D >= 2 ? y2glob - mins[1] : y2glob;
  z1 = D >= 3 ? z1glob - mins[2] : z1glob;
  z2 = D >= 3 ? z2glob - mins[2] : z2glob;

  i1  = D >= 1 ? static_cast<int>(floor( x1 ) ) : 0;
  i2  = D >= 1 ? static_cast<int>(floor( x2 ) ) : 0;
  j1  = D >= 2 ? static_cast<int>(floor( y1 ) ) : 0;
  j2  = D >= 2 ? static_cast<int>(floor( y2 ) ) : 0;
  k1  = D >= 3 ? static_cast<int>(floor( z1 ) ) : 0;
  k2  = D >= 3 ? static_cast<int>(floor( z2 ) ) : 0;


  // index checking
  if(D >= 1 ) assert(i1 >= -3 && i1 < static_cast<int>(tile.mesh_lengths[0]+3)) ;
  if(D >= 2 ) assert(j1 >= -3 && j1 < static_cast<int>(tile.mesh_lengths[1]+3)) ;
  if(D >= 3 ) assert(k1 >= -3 && k1 < static_cast<int>(tile.mesh_lengths[2]+3)) ;

  if (D >= 1) assert(i2 >= -3 && i2 < static_cast<int>(tile.mesh_lengths[0]+3));
  if (D >= 2) assert(j2 >= -3 && j2 < static_cast<int>(tile.mesh_lengths[1]+3));
  if (D >= 3) assert(k2 >= -3 && k2 < static_cast<int>(tile.mesh_lengths[2]+3));


  // relay point; +1 is equal to +\Delta x
  xr = min( (real_long)(min(i1,i2)+1), max( (real_long)max(i1,i2), static_cast<real_long>(0.5)*(x1+x2) ) );
  yr = min( (real_long)(min(j1,j2)+1), max( (real_long)max(j1,j2), static_cast<real_long>(0.5)*(y1+y2) ) );
  zr = min( (real_long)(min(k1,k2)+1), max( (real_long)max(k1,k2), static_cast<real_long>(0.5)*(z1+z2) ) );


  // -q to include -j in the Ampere's equation
  Fx1 = +q*(xr - x1);
  Fy1 = +q*(yr - y1);
  Fz1 = +q*(zr - z1);
  
  Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
  Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
  Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

  Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
  Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;

  Fx2 = +q*(x2-xr);
  Fy2 = +q*(y2-yr);
  Fz2 = +q*(z2-zr);

  // jx
  if (D >= 1) yee.jx(i1,  j1,   k1)   += Fx1 * (1.0-Wy1) * (1.0-Wz1);
  if (D >= 2) yee.jx(i1,  j1+1, k1)   += Fx1 * Wy1       * (1.0-Wz1);
  if (D >= 3) yee.jx(i1,  j1,   k1+1) += Fx1 * (1.0-Wy1) * Wz1;
  if (D >= 3) yee.jx(i1,  j1+1, k1+1) += Fx1 * Wy1     * Wz1;

  if (D >= 1) yee.jx(i2,  j2,   k2)   += Fx2 * (1.0-Wy2) * (1.0-Wz2);
  if (D >= 2) yee.jx(i2,  j2+1, k2)   += Fx2 * Wy2       * (1.0-Wz2);
  if (D >= 3) yee.jx(i2,  j2,   k2+1) += Fx2 * (1.0-Wy2) * Wz2;
  if (D >= 3) yee.jx(i2,  j2+1, k2+1) += Fx2 * Wy2     * Wz2;

  // jy
  if (D >= 1) yee.jy(i1,  j1,   k1)   += Fy1 * (1.0-Wx1) * (1.0-Wz1);
  if (D >= 2) yee.jy(i1+1,j1,   k1)   += Fy1 * Wx1       * (1.0-Wz1);
  if (D >= 3) yee.jy(i1  ,j1,   k1+1) += Fy1 * (1.0-Wx1) * Wz1;
  if (D >= 3) yee.jy(i1+1,j1,   k1+1) += Fy1 * Wx1     * Wz1;
  
  if (D >= 1) yee.jy(i2,  j2,   k2)   += Fy2 * (1.0-Wx2) * (1.0-Wz2);
  if (D >= 2) yee.jy(i2+1,j2,   k2)   += Fy2 * Wx2       * (1.0-Wz2);
  if (D >= 3) yee.jy(i2,  j2,   k2+1) += Fy2 * (1.0-Wx2) * Wz2;
  if (D >= 3) yee.jy(i2+1,j2,   k2+1) += Fy2 * Wx2     * Wz2;
                        
  // jz
  yee.jz(i1,  j1,   k1)   += Fz1 * (1.0-Wx1) * (1.0-Wy1);
  yee.jz(i1+1,j1,   k1)   += Fz1 * Wx1       * (1.0-Wy1);
  yee.jz(i1,  j1+1, k1)   += Fz1 * (1.0-Wx1) * Wy1;
  yee.jz(i1+1,j1+1, k1)   += Fz1 * Wx1       * Wy1;

  yee.jz(i2,  j2,   k2)   += Fz2 * (1.0-Wx2) * (1.0-Wy2);
  yee.jz(i2+1,j2,   k2)   += Fz2 * Wx2       * (1.0-Wy2);
  yee.jz(i2,  j2+1, k2)   += Fz2 * (1.0-Wx2) * Wy2;
  yee.jz(i2+1,j2+1, k2)   += Fz2 * Wx2       * Wy2;

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
    int nparts = container.size();

    // initialize pointers to particle arrays
    real_prtcl* loc[3];
    for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    real_prtcl* vel[3];
    for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    real_long c = tile.cfl;
    real_long loc0n, loc1n, loc2n, vel0n, vel1n, vel2n;
    real_long x0, y0, z0, gamma, tfrac, walloc0;
    real_long xcolis, ycolis, zcolis;

    for(int n=n1; n<n2; n++) {

      // left side of the wall boundary
      if( static_cast<real_long>(loc[0][n]) < walloc) {

        loc0n = static_cast<real_long>( loc[0][n] );
        loc1n = static_cast<real_long>( loc[1][n] );
        loc2n = static_cast<real_long>( loc[2][n] );
                                                     
        vel0n = static_cast<real_long>( vel[0][n] );
        vel1n = static_cast<real_long>( vel[1][n] );
        vel2n = static_cast<real_long>( vel[2][n] );


        // unwind wall location
        walloc0 = walloc - betawall*c;

        gamma = sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);
        x0 = loc0n - vel0n/gamma*c;
        y0 = loc1n;
        z0 = loc2n;

        // this can not be reflected particle; remove 
        // equals to right boundary outflow as they wrap particles here
        if(walloc0 - x0 > c) {
          //std::cout << "reflected?" << n << "x:" << loc[0][n] << "\n";
          to_be_deleted.push_back(n);
          continue;
        }

        // compute crossing point
        //tfrac = std::min(
        //    std::abs((x0-walloc0)/(betawall*c - vel0n/gamma*c)), 
        //    (real_long)1.0
        //    );
        tfrac = std::abs((x0-walloc0)/(betawall*c - vel0n/gamma*c));
        if (tfrac > 1.0) {
           // std::cout << "current up to intersection point:\n"
           //   << " x0      " << x0
           //   << " vel0n   " << vel0n/gamma
           //   << " tfrac   " << tfrac
           //   << " walloc  " << walloc
           //   << " walloc0 " << walloc0
           //   << "\n";

            to_be_deleted.push_back(n);
            continue;
        }

        xcolis = x0 + vel0n/gamma*c*tfrac;
        ycolis = y0;
        zcolis = z0;

        /*
        std::cout << "current up to intersection point:\n"
          << " x0      " << x0
          << " xcolis  " << xcolis
          << " tfrac   " << tfrac
          << " walloc  " << walloc
          << " walloc0 " << walloc0
          << "\n";
        */
          
        // deposit current up to intersection point
        zigzag(tile, xcolis, ycolis, zcolis, x0, y0, z0, container.q);

        // reset particle momentum, getting kick from the wall
        vel0n = gammawall*gammawall*gamma
          *(2.*betawall - vel0n/gamma*(1.0+betawall*betawall));

        // recompute gamma
        gamma = sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

        tfrac = std::min(
            std::abs((vel0n-xcolis)/std::max(std::abs(vel0n-x0), (real_long)1.0e-6)), 
            (real_long)1.0
            );

        // move particle from the location of the wall with new velocity
        loc0n = xcolis + vel0n/gamma*c * tfrac;
        loc1n = ycolis;
        loc2n = zcolis;

        /*
        std::cout << "current cleaning behind wall:\n"
          << " xcolis  " << xcolis
          << " x0      " << loc0n - vel0n/gamma*c
          << " tfrac   " << tfrac
          << " gamma   " << gamma
          << "\n";
        */

        // clean up the part of trajectory behind the wall
        // that will be added by the deposition routine that unwinds
        // the particle location by a full time step.
        zigzag(
            tile, 
            xcolis, ycolis, zcolis,
            loc0n - vel0n/gamma*c,
            loc1n - vel1n/gamma*c,
            loc2n - vel2n/gamma*c,
            -container.q);

        // lastly; store particle back to the container
        loc[0][n] = static_cast<real_prtcl>( loc0n );
        loc[1][n] = static_cast<real_prtcl>( loc1n );
        loc[2][n] = static_cast<real_prtcl>( loc2n );

        vel[0][n] = static_cast<real_prtcl>( vel0n );
        vel[1][n] = static_cast<real_prtcl>( vel1n );
        vel[2][n] = static_cast<real_prtcl>( vel2n );
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
  if(walloc < maxs[0]) {
    auto& yee = tile.get_yee();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0]);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) {
      for(int i=-3; i<=iw; i++) {

        // transverse components of electric field to zero (only parallel comp allowed)
        yee.ey(i,j,k) = 0.0;
        yee.ez(i,j,k) = 0.0;

        // parallel comp to zero (conductor repels B field)
        yee.bx(i,j,k) = 0.0;

        // clean all current behind piston head
        yee.jx(i,j,k) = 0.0;
        yee.jy(i,j,k) = 0.0;
        yee.jz(i,j,k) = 0.0;
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
    auto& yee = tile.get_yee();

    // wall location 
    auto iw = static_cast<int>(walloc - mins[0]);
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero to make this conductor
    for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
    for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
    for(int i=-3; i<=iw; i++) {

      // transverse components of electric field to zero (only parallel comp allowed)
      yee.ey(i,j,k) = 0.0;
      yee.ez(i,j,k) = 0.0;

      // parallel comp to zero (conductor repels B field)
      yee.bx(i,j,k) = 0.0;

      // clean all current behind piston head
      yee.jx(i,j,k) = 0.0;
      yee.jy(i,j,k) = 0.0;
      yee.jz(i,j,k) = 0.0;
    }

  } // end if inside piston

}


//--------------------------------------------------
// explicit template instantiation

//template class pic::Piston<1>; // 1D3V
template class pic::Piston<2>; // 2D3V
template class pic::Piston<3>; // 3D3V

