#include "boris_drag.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;


template<size_t D, size_t V>
double pic::BorisPusherDrag<D,V>::kn(double x)
{
  if (temp == 0.0) return 1.0; // Thompson limit

  double sig;
  sig = (1.0 - 4.0/x - 8.0/x/x)*log(1.+x) + 0.5 + 8.0/x - 1.0/(2.0*(1. + x)*(1. + x));
  return (3.0/4.0)*sig/x;
}



template<size_t D, size_t V>
void pic::BorisPusherDrag<D,V>::push_container(
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
  real_long uxt, uyt, uzt;
  real_long u1, v1, w1;
  real_long g, f, ginv, kncorr, gamt, ut;
  real_long dragx, dragy, dragz, dragv;
  real_long thr;

  // maximum drag force experienced by particle
  real_long dragthr = 0.1; 

  real_long c = cfl;
  real_long cinv = 1.0/c;

  // charge (sign only)
  real_long qm = sign(container.q);

  real_long vel0n, vel1n, vel2n;

  for(int n=n1; n<n2; n++) {

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


    //--------------------------------------------------
    // addition of drag (gamma at half time step)

    // u at t + dt/2
    uxt  = (u0*cinv + vel0n)*0.5;
    uyt  = (v0*cinv + vel1n)*0.5;
    uzt  = (w0*cinv + vel2n)*0.5;
    ut   = sqrt(uxt*uxt + uyt*uyt + uzt*uzt);
    gamt = sqrt(1.0 + ut*ut);

    // subtract drag with Klein-Nishina reduction
    // A g^2 beta = A g^2 u/g = A g u
    kncorr = kn(3.0*gamt*temp);

    // drag components
    dragx = c*drag*kncorr*gamt*gamt*(uxt/gamt);
    dragy = c*drag*kncorr*gamt*gamt*(uyt/gamt);
    dragz = c*drag*kncorr*gamt*gamt*(uzt/gamt);

    // limit drag to maximum of dragthr of velocity
    dragv = sqrt(dragx*dragx + dragy*dragy + dragz*dragz)/ut;
    thr = 1.0;
    if (dragv > dragthr) thr = dragthr/dragv;

    // apply drag
    vel[0][n] = static_cast<real_prtcl>( u0*cinv - thr*dragx );
    vel[1][n] = static_cast<real_prtcl>( v0*cinv - thr*dragy );
    vel[2][n] = static_cast<real_prtcl>( w0*cinv - thr*dragz );


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

template class pic::BorisPusherDrag<1,3>; // 1D3V
template class pic::BorisPusherDrag<2,3>; // 2D3V
template class pic::BorisPusherDrag<3,3>; // 3D3V

