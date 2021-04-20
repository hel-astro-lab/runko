#include "piston_z.h"

#include <cmath> 
#include <cassert>
#include <algorithm>

#include "../../tools/signum.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


using std::min;
using std::max;



template<size_t D>
void pic::PistonZdir<D>::zigzag(
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
  auto& yee = tile.get_yee();

  const auto mins = tile.mins;
  const double c = tile.cfl;    // speed of light

  // normalized location w.r.t. tile
  float_m x1 = D >= 1 ? x1glob - mins[0] : x1glob;
  float_m x2 = D >= 1 ? x2glob - mins[0] : x2glob;

  float_m y1 = D >= 2 ? y1glob - mins[1] : y1glob;
  float_m y2 = D >= 2 ? y2glob - mins[1] : y2glob;

  float_m z1 = D >= 3 ? z1glob - mins[2] : z1glob;
  float_m z2 = D >= 3 ? z2glob - mins[2] : z2glob;

  int i1  = D >= 1 ? floor( x1 ) : 0;
  int i2  = D >= 1 ? floor( x2 ) : 0;

  int j1  = D >= 2 ? floor( y1 ) : 0;
  int j2  = D >= 2 ? floor( y2 ) : 0;

  int k1  = D >= 3 ? floor( z1 ) : 0;
  int k2  = D >= 3 ? floor( z2 ) : 0;

  // relay point; +1 is equal to +\Delta x
  float_m xr = min( float_m(min(i1,i2)+1), max( float_m(max(i1,i2)), float_m(0.5*(x1+x2)) ) );
  float_m yr = min( float_m(min(j1,j2)+1), max( float_m(max(j1,j2)), float_m(0.5*(y1+y2)) ) );
  float_m zr = min( float_m(min(k1,k2)+1), max( float_m(max(k1,k2)), float_m(0.5*(z1+z2)) ) );

  //--------------------------------------------------
  // +q since - sign is already included in the Ampere's equation
  float_m Fx1 = +q*(xr - x1);
  float_m Fy1 = +q*(yr - y1);
  float_m Fz1 = +q*(zr - z1);
  
  float_m Fx2 = +q*(x2 - xr);
  float_m Fy2 = +q*(y2 - yr);
  float_m Fz2 = +q*(z2 - zr);

  float_m Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
  float_m Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
  float_m Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

  float_m Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
  float_m Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  float_m Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;


  //--------------------------------------------------
  // jx
  if(D>=1) atomic_add( yee.jx(i1  , j1  , k1  ), Fx1*(1.0f-Wy1)*(1.0f-Wz1) );
  if(D>=2) atomic_add( yee.jx(i1  , j1+1, k1  ), Fx1*Wy1       *(1.0f-Wz1) );
  if(D>=3) atomic_add( yee.jx(i1  , j1  , k1+1), Fx1*(1.0f-Wy1)*Wz1        );
  if(D>=3) atomic_add( yee.jx(i1  , j1+1, k1+1), Fx1*Wy1       *Wz1        );

  if(D>=1) atomic_add( yee.jx(i2  , j2  , k2  ), Fx2*(1.0f-Wy2)*(1.0f-Wz2) );
  if(D>=2) atomic_add( yee.jx(i2  , j2+1, k2  ), Fx2*Wy2       *(1.0f-Wz2) );
  if(D>=3) atomic_add( yee.jx(i2  , j2  , k2+1), Fx2*(1.0f-Wy2)*Wz2        );
  if(D>=3) atomic_add( yee.jx(i2  , j2+1, k2+1), Fx2*Wy2       *Wz2        );

  // jy
  if(D>=1) atomic_add( yee.jy(i1  , j1  , k1  ), Fy1*(1.0f-Wx1)*(1.0f-Wz1) );
  if(D>=2) atomic_add( yee.jy(i1+1, j1  , k1  ), Fy1*Wx1       *(1.0f-Wz1) );
  if(D>=3) atomic_add( yee.jy(i1  , j1  , k1+1), Fy1*(1.0f-Wx1)*Wz1        );
  if(D>=3) atomic_add( yee.jy(i1+1, j1  , k1+1), Fy1*Wx1       *Wz1        );

  if(D>=1) atomic_add( yee.jy(i2  , j2  , k2  ), Fy2*(1.0f-Wx2)*(1.0f-Wz2) );
  if(D>=2) atomic_add( yee.jy(i2+1, j2  , k2  ), Fy2*Wx2       *(1.0f-Wz2) );
  if(D>=3) atomic_add( yee.jy(i2  , j2  , k2+1), Fy2*(1.0f-Wx2)*Wz2        );
  if(D>=3) atomic_add( yee.jy(i2+1, j2  , k2+1), Fy2*Wx2       *Wz2        );

  // jz
  if(D>=1) atomic_add( yee.jz(i1  , j1  , k1  ), Fz1*(1.0f-Wx1)*(1.0f-Wy1) );
  if(D>=2) atomic_add( yee.jz(i1+1, j1  , k1  ), Fz1*Wx1       *(1.0f-Wy1) );
  if(D>=3) atomic_add( yee.jz(i1  , j1+1, k1  ), Fz1*(1.0f-Wx1)*Wy1        );
  if(D>=3) atomic_add( yee.jz(i1+1, j1+1, k1  ), Fz1*Wx1       *Wy1        );

  if(D>=1) atomic_add( yee.jz(i2  , j2  , k2  ), Fz2*(1.0f-Wx2)*(1.0f-Wy2) );
  if(D>=1) atomic_add( yee.jz(i2+1, j2  , k2  ), Fz2*Wx2       *(1.0f-Wy2) );
  if(D>=2) atomic_add( yee.jz(i2  , j2+1, k2  ), Fz2*(1.0f-Wx2)*Wy2        );
  if(D>=2) atomic_add( yee.jz(i2+1, j2+1, k2  ), Fz2*Wx2       *Wy2        );

}


template<size_t D>
void pic::PistonZdir<D>::solve( pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const auto mins = tile.mins;
  const auto maxs = tile.maxs;

  // unwind wall location
  double wz0 = wallocz; // - wdirz*betawallz*c;
  bool tile_between_z = (mins[2] <= wz0 && wz0 <= maxs[2]);


  // if wall crossing is not inside tile, skip further analysis
  if( !tile_between_z ) return;

  for(auto&& con : tile.containers) {
    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge

    UniIter::iterate([=] DEVCALLABLE (
                size_t n, 
                pic::Tile<D>& tile,
                pic::ParticleContainer<D>& con
                ){

      //particle location and velocity
      double x1 = con.loc(0,n);
      double y1 = con.loc(1,n);
      double z1 = con.loc(2,n);

      bool do_reflect = false;
      if( wdir > 0.0 && z1 < wallocz ) do_reflect = true;
      if( wdir < 0.0 && z1 > wallocz ) do_reflect = true;

      // skip if we dont reflect
      if(do_reflect) {

        //--------------------------------------------------
        double u1 = con.vel(0,n);
        double v1 = con.vel(1,n);
        double w1 = con.vel(2,n);
        double gam = sqrt(1.0 + u1*u1 + v1*v1 + w1*w1);

        // unwind prtcl location; 
        // NOTE: this could be approximated by only treating the refl. component
        double x0 = x1 - c*u1/gam;
        double y0 = y1 - c*v1/gam;
        double z0 = z1 - c*w1/gam;

        //--------------------------------------------------
        // compute crossing point

        // time fraction between previos position and wall
        double dt = std::min(1.0, std::abs((z0-wallocz)/(-wdir*c*w1/gam + EPS)));
        double xcol = x0 + c*dt*u1/gam;
        double ycol = y0 + c*dt*v1/gam;
        double zcol = z0 + c*dt*w1/gam;

        // deposit current upto intersection w/ wall
        zigzag(tile, x0, y0, z0, xcol, ycol, zcol, q);

        //--------------------------------------------------
        // perform reflection

        // reflect particle getting kick from the wall
        //w1 = gamw*gamw*gam*(2.*wdirz*betawallz - w/gam*(1.0 + betawallz*betawallz));
        w1 = - w1;

        // time fraction between current position and collision point
        dt = std::min(1.0, std::abs((z1-zcol)/(z1 - z0 + EPS)) );
        // NOTE: this could be approximated by only treating the refl. component
        double xnew = xcol + c*dt*u1/gam;
        double ynew = ycol + c*dt*v1/gam;
        double znew = zcol + c*dt*w1/gam;

        con.loc(0,n) = xnew;
        con.loc(1,n) = ynew;
        con.loc(2,n) = znew;

        //con.vel(0,n) = u1;
        //con.vel(1,n) = v1;
        con.vel(2,n) = w1;

        // clean up the part of trajectory behind the wall
        // that will be added by the deposition routine that unwinds
        // the particle location by a full time step.
        zigzag(tile, 
            xnew - c*u1/gam,
            ynew - c*v1/gam,
            znew - c*w1/gam,
            xcol, ycol, zcol, 
            -q);

        }


      }, 
        con.size(), 
        tile,
        con);

    UniIter::sync();
  }


#ifdef GPU
  nvtxRangePop();
#endif
}


template<>
void pic::PistonZdir<3>::field_bc(
    pic::Tile<3>& tile)
{

  // skip if piston head is not inside tile boundaries
  const auto mins = tile.mins;
  const auto maxs = tile.maxs;
  auto& yee = tile.get_yee();

  // walloc to relative tile units
  int kw = wallocz - mins[2];

  // make left side of piston conductor
  if(wdir > 0 && wallocz < maxs[2]) {

    // limit wall location 
    if(kw > tile.mesh_lengths[2]) kw = tile.mesh_lengths[2];

    // set transverse directions to zero to make this conductor
    for(int k=-3; k<=kw; k++) 
    for(int j=-3; j<tile.mesh_lengths[1]+3; j++) 
    for(int i=-3; i<tile.mesh_lengths[0]+3; i++) {

      // transverse components of electric field to zero (only parallel comp allowed)
      yee.ey(i,j,k) = 0.0;
      yee.ez(i,j,k) = 0.0;

      // clean all current behind piston head
      //yee.jx(i,j,k) = 0.0;
      //yee.jy(i,j,k) = 0.0;
      //yee.jz(i,j,k) = 0.0;
    }
  } else if(wdir < 0 && wallocz > mins[2]) {

    // limit wall location
    if(kw < -3) kw = -3;

    // set transverse directions to zero to make this conductor
    for(int k=kw; k<tile.mesh_lengths[2]+3; k++) 
    for(int j=-3; j<tile.mesh_lengths[1]+3; j++) 
    for(int i=-3; i<tile.mesh_lengths[0]+3; i++) {

      // transverse components of electric field to zero (only parallel comp allowed)
      yee.ey(i,j,k) = 0.0;
      yee.ez(i,j,k) = 0.0;

      // clean all current behind piston head
      //yee.jx(i,j,k) = 0.0;
      //yee.jy(i,j,k) = 0.0;
      //yee.jz(i,j,k) = 0.0;
    }
  }
}


template<size_t D>
void pic::PistonZdir<D>::clean_prtcls( pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const auto mins = tile.mins;
  const auto maxs = tile.maxs;

  // unwind wall location
  double wz0 = wallocz; // - wdirz*betawallz*c;
  bool tile_between_z = (mins[2] <= wz0 && wz0 <= maxs[2]);

  // if wall crossing is not inside tile, skip further analysis
  if( !tile_between_z ) return;

  // outflowing particles
  std::vector<int> to_be_deleted;


  for(auto&& con : tile.containers) {
    const double c = tile.cfl;    // speed of light
    to_be_deleted.clear();

    for(int n=0; n<con.size(); n++) {

      //particle location and velocity
      double z1 = con.loc(2,n);

      bool behind_wall = false;
      if( wdir > 0.0 && z1 < wallocz ) behind_wall = true;
      if( wdir < 0.0 && z1 > wallocz ) behind_wall = true;

      // skip if we dont reflect
      if(behind_wall) to_be_deleted.push_back(n);
    }

    con.delete_particles(to_be_deleted);
  }

}

//template class pic::Piston<1>; // 1D3V
//template class pic::PistonZdir<2>; // 2D3V
template class pic::PistonZdir<3>; // 3D3V
