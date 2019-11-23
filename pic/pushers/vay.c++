#include "vay.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::VayPusher<D,V>::push_container(
    pic::ParticleContainer& container, 
    double cfl) 
{
  int nparts = container.size();

  // initialize pointers to particle arrays
  float_tp* loc[3];
  for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

  float_tp* vel[3];
  for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );


  double_t ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
  double_t bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

  // make sure E and B tmp arrays are of correct size
  if(container.Epart.size() != (size_t)3*nparts)
    container.Epart.resize(3*nparts);
  if(container.Bpart.size() != (size_t)3*nparts)
    container.Bpart.resize(3*nparts);

  float_tp *ex, *ey, *ez, *bx, *by, *bz;
  ex = &( container.Epart[0*nparts] );
  ey = &( container.Epart[1*nparts] );
  ez = &( container.Epart[2*nparts] );

  bx = &( container.Bpart[0*nparts] );
  by = &( container.Bpart[1*nparts] );
  bz = &( container.Bpart[2*nparts] );

  // loop over particles
  int n1 = 0;
  int n2 = nparts;

  double_t u0, v0, w0;
  double_t u1, v1, w1;
  double_t g, f;
  double_t ustar, sig, tx, ty, tz, vx0, vy0, vz0;

  double_t c = cfl;
  double_t cinv = 1.0/c;

  // charge (sign only)
  double_t qm = sign(container.q);

  double_t vel0n, vel1n, vel2n;

  // add division by m_s to simulate multiple species

  //TODO: SIMD
  for(int n=n1; n<n2; n++) {
    vel0n = static_cast<double_t>( vel[0][n] );
    vel1n = static_cast<double_t>( vel[1][n] );
    vel2n = static_cast<double_t>( vel[2][n] );

    // read particle-specific fields
    ex0 = static_cast<double_t>( ex[n]*(0.5*qm) );
    ey0 = static_cast<double_t>( ey[n]*(0.5*qm) );
    ez0 = static_cast<double_t>( ez[n]*(0.5*qm) );

    bx0 = static_cast<double_t>( bx[n]*(0.5*qm*cinv) );
    by0 = static_cast<double_t>( by[n]*(0.5*qm*cinv) );
    bz0 = static_cast<double_t>( bz[n]*(0.5*qm*cinv) );

    //--------------------------------------------------
    // Vay algorithm
      
    // gamma^-1
    g = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

    vx0 = c*vel0n*g;
    vy0 = c*vel1n*g;
    vz0 = c*vel2n*g;

    // u' (cinv is already multiplied into B)
    u1 = c*vel0n + 2.0*ex0 + vy0*bz0 - vz0*by0;
    v1 = c*vel1n + 2.0*ey0 + vz0*bx0 - vx0*bz0;
    w1 = c*vel2n + 2.0*ez0 + vx0*by0 - vy0*bx0;
    
    // gamma(u')
    ustar = cinv*(u1*bx0+v1*by0+w1*bz0);
		sig = cinv*cinv*( c*c + u1*u1+ v1*v1+ w1*w1) - (bx0*bx0+ by0*by0+ bz0*bz0);
		g = 1.0/sqrt( 0.5*(sig + sqrt(sig*sig + 4.0*(bx0*bx0 + by0*by0 + bz0*bz0 + ustar*ustar))));

    tx = bx0*g;
	  ty = by0*g;
	  tz = bz0*g;
	  f = 1.0/(1.0+ tx*tx+ ty*ty+ tz*tz);
	  
    // final step 
		u0 = f*(u1 + (u1*tx + v1*ty + w1*tz)*tx + v1*tz - w1*ty);
		v0 = f*(v1 + (u1*tx + v1*ty + w1*tz)*ty + w1*tx - u1*tz);
		w0 = f*(w1 + (u1*tx + v1*ty + w1*tz)*tz + u1*ty - v1*tx);

    // normalized 4-velocity advance
    vel[0][n] = static_cast<float_tp>( u0*cinv );
    vel[1][n] = static_cast<float_tp>( v0*cinv );
    vel[2][n] = static_cast<float_tp>( w0*cinv );

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    g = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n]*g*c;
  }
}



//--------------------------------------------------
// explicit template instantiation

template class pic::VayPusher<1,3>; // 1D3V
template class pic::VayPusher<2,3>; // 2D3V
template class pic::VayPusher<3,3>; // 3D3V

