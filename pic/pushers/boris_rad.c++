#include "boris_rad.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::BorisPusherRad<D,V>::push_container(
    pic::ParticleContainer<D>& container, 
    double cfl)
{

  int nparts = container.size();

  // initialize pointers to particle arrays
  real_prtcl* loc[3];
  for( int i=0; i<3; i++)
    loc[i] = &( container.loc(i,0) );

  real_prtcl* vel[3];
  for( int i=0; i<3; i++)
    vel[i] = &( container.vel(i,0) );


  real_long ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
  real_long bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

  // make sure E and B tmp arrays are of correct size
  if(container.Epart.size() != (size_t)3*nparts)
    container.Epart.resize(3*nparts);
  if(container.Bpart.size() != (size_t)3*nparts)
    container.Bpart.resize(3*nparts);

  real_prtcl *ex, *ey, *ez, *bx, *by, *bz;
  ex = &( container.Epart[0*nparts] );
  ey = &( container.Epart[1*nparts] );
  ez = &( container.Epart[2*nparts] );

  bx = &( container.Bpart[0*nparts] );
  by = &( container.Bpart[1*nparts] );
  bz = &( container.Bpart[2*nparts] );

  // loop over particles
  int n1 = 0;
  int n2 = nparts;

  real_long u0, v0, w0;
  real_long uxt;
  real_long u1, v1, w1;
  real_long g, f, ginv;
  real_long pressx, pressy, pressz;
  real_long gamx;

  real_long c = cfl;
  real_long cinv = 1.0/c;

  // charge (sign only)
  real_long qm = sign(container.q);

  real_long loc0n;
  real_long vel0n, vel1n, vel2n;

  for(int n=n1; n<n2; n++) {

    loc0n = static_cast<real_long>( loc[0][n] ); // x position

    vel0n = static_cast<real_long>( vel[0][n] );
    vel1n = static_cast<real_long>( vel[1][n] );
    vel2n = static_cast<real_long>( vel[2][n] );

    //--------------------------------------------------
    // Boris algorithm

    // read particle-specific fields
    ex0 = static_cast<real_long>( ex[n] )*0.5*qm;
    ey0 = static_cast<real_long>( ey[n] )*0.5*qm;
    ez0 = static_cast<real_long>( ez[n] )*0.5*qm;

    bx0 = static_cast<real_long>( bx[n] )*0.5*qm*cinv;
    by0 = static_cast<real_long>( by[n] )*0.5*qm*cinv;
    bz0 = static_cast<real_long>( bz[n] )*0.5*qm*cinv;

    // first half electric acceleration
    u0 = c*vel0n + ex0;
    v0 = c*vel1n + ey0;
    w0 = c*vel2n + ez0;

    // first half magnetic rotation
    g = c/sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    bx0 *= g;
    by0 *= g;
    bz0 *= g;

    f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
    u1 = (u0 + v0*bz0 - w0*by0)*f;
    v1 = (v0 + w0*bx0 - u0*bz0)*f;
    w1 = (w0 + u0*by0 - v0*bx0)*f;

    // second half of magnetic rotation & electric acceleration
    u0 = u0 + v1*bz0 - w1*by0 + ex0;
    v0 = v0 + w1*bx0 - u1*bz0 + ey0;
    w0 = w0 + u1*by0 - v1*bx0 + ez0;


    // rad pressure components; apply pressure only to particles behind the rad beam front
    if(loc0n <= beam_locx){

        // addition of radiation pressure (gamma at half time step)
        // u at t + dt/2
        uxt  = (u0*cinv + vel0n)*0.5;
        // uyt  = (v0*cinv + vel1n)*0.5;
        // uzt  = (w0*cinv + vel2n)*0.5;
        // ut   = sqrt(uxt*uxt + uyt*uyt + uzt*uzt);


        // betax component of the prtcl velocity
        // gamt = sqrt(1.0 + ut*ut);
        // betax = uxt/gamt;
        //gamx  = 1.0/sqrt(1.0 - betax*betax):
        gamx = sqrt(1.0 + uxt*uxt);

        pressx = c*drag/gamx/gamx;
        pressy = 0;
        pressz = 0;
    } else {
        pressx = 0;
        pressy = 0;
        pressz = 0;
    }


    // apply pressure
    vel[0][n] = static_cast<real_prtcl>( u0*cinv + pressx );
    vel[1][n] = static_cast<real_prtcl>( v0*cinv + pressy );
    vel[2][n] = static_cast<real_prtcl>( w0*cinv + pressz );

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = 1.0/sqrt(1.0 + 
        vel[0][n]*vel[0][n] +
        vel[1][n]*vel[1][n] +
        vel[2][n]*vel[2][n]);

    for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n]*ginv*c;
  }
}


//template<size_t D, size_t V>
//void pic::BorisPusherDrag<D,V>::solve(
//    pic::Tile<D>& tile)
//{
//
//  for(auto&& container : tile.containers)
//    push_container(container, tile.cfl);
//
//}



//--------------------------------------------------
// explicit template instantiation

template class pic::BorisPusherRad<1,3>; // 1D3V
template class pic::BorisPusherRad<2,3>; // 2D3V
template class pic::BorisPusherRad<3,3>; // 3D3V
