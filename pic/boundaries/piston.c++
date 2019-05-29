#include "piston.h"
#include "../../tools/signum.h"

#include <cmath> 
#include <cassert>
#include <algorithm>

using std::min;
using std::max;

//TODO: FIXME directy copy from depositer/zigzag.c++
// use friending of the class or inherit this?
template<size_t D>
void pic::Piston<D>::zigzag(
    pic::Tile<D>& tile,
    double x2glob, 
    double y2glob, 
    double z2glob, 
    double x1glob, 
    double y1glob, 
    double z1glob, 
    double q)
{
  int i1,i2,j1,j2,k1,k2;

  double x1,x2,y1,y2,z1,z2;
  double xr,yr,zr;

  double Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
  double Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;

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
  xr = min( (double)(min(i1,i2)+1), max( (double)max(i1,i2), static_cast<double>(0.5)*(x1+x2) ) );
  yr = min( (double)(min(j1,j2)+1), max( (double)max(j1,j2), static_cast<double>(0.5)*(y1+y2) ) );
  zr = min( (double)(min(k1,k2)+1), max( (double)max(k1,k2), static_cast<double>(0.5)*(z1+z2) ) );


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
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    double c = tile.cfl;
    //double qm = sign(container.q);

    double x0, y0, z0, gamma, tfrac, walloc0;
    double xcolis, ycolis, zcolis;
    for(int n=n1; n<n2; n++) {

      // left side of the wall boundary
      if(loc[0][n] < walloc) {

        // this can not be reflected particle; remove 
        // equals to right boundary outflow as they wrap particles here
        if(walloc - loc[0][n] > 1.0) {
          //std::cout << "reflected?" << n << "x:" << loc[0][n] << "\n";
          to_be_deleted.push_back(n);
          continue;
        }

        gamma = sqrt(1.0
            + vel[0][n]*vel[0][n]
            + vel[1][n]*vel[1][n]
            + vel[2][n]*vel[2][n]);

        x0 = loc[0][n] - vel[0][n]/gamma*c;
        y0 = loc[1][n];
        z0 = loc[2][n];

        // unwind wall location
        walloc0 = walloc - betawall*c;

        // compute crossing point
        tfrac = abs((x0-walloc0)/(betawall*c - vel[0][n]/gamma*c));
        xcolis = x0 + vel[0][n]/gamma*c*tfrac;
        ycolis = y0;
        zcolis = z0;
          
        // deposit current up to intersection point
        zigzag(tile, xcolis, ycolis, zcolis, x0, y0, z0, container.q);

        // reset particle momentum, getting kick from the wall
        vel[0][n] = gammawall*gammawall*gamma
          *(2.*betawall - vel[0][n]/gamma*(1.0+betawall*betawall));

        tfrac = std::min( std::abs((vel[0][n]-xcolis)/std::max(abs(vel[0][n]-x0), 1.0e-9)), 1.0);

        // move particle from the location of the wall with new velocity
        loc[0][n] = xcolis + vel[0][n]/gamma*c * tfrac;
        loc[1][n] = ycolis;
        loc[2][n] = zcolis;

        // clean up the part of trajectory behind the wall
        // that will be added by the deposition routine that unwinds
        // the particle location by a full time step.
        zigzag(
            tile, 
            xcolis, ycolis, zcolis,
            loc[0][n] - vel[0][n]/gamma*c,
            loc[1][n] - vel[1][n]/gamma*c,
            loc[2][n] - vel[2][n]/gamma*c,
            -container.q);
      }
    }

    // process outflown particles
    container.delete_particles(to_be_deleted);

  } // end of loop over species

  return;
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
    int iw = walloc - mins[0]; 
    if(iw > static_cast<int>(tile.mesh_lengths[0])) iw = tile.mesh_lengths[0];

    // set transverse directions to zero
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<=iw; i++) {
        yee.ey(i,j,k) = 0.0;
        yee.ez(i,j,k) = 0.0;
      }
    }
  }

}


//--------------------------------------------------
// explicit template instantiation

//template class pic::Piston<1>; // 1D3V
template class pic::Piston<2>; // 2D3V
//template class pic::Piston<3>; // 3D3V

