#include "higuera_cary.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::HigueraCaryPusher<D,V>::push_container(
    pic::ParticleContainer<D>& container, 
    double cfl) 
{
  int nparts = container.size();

  // initialize pointers to particle arrays
  real_prtcl* loc[3];
  for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

  real_prtcl* vel[3];
  for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );


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

  real_long c = cfl;
  real_long cinv = 1.0/c;

  // half charge-to-mass ratio (sign only because fields are in units of q)
  real_long qm2 = 0.5*sign(container.q)/container.m;

  real_long vel0n, vel1n, vel2n;
  real_long u0, v0, w0;
  real_long u1, v1, w1;
  real_long g2, b2;
  real_long ginv, f;


  for(int n=n1; n<n2; n++) {
    vel0n = static_cast<real_long>( vel[0][n] );
    vel1n = static_cast<real_long>( vel[1][n] );
    vel2n = static_cast<real_long>( vel[2][n] );

    // read particle-specific fields
    ex0 = static_cast<real_long>( (ex[n] + this->get_ex_ext(0,0,0))*qm2 );
    ey0 = static_cast<real_long>( (ey[n] + this->get_ey_ext(0,0,0))*qm2 );
    ez0 = static_cast<real_long>( (ez[n] + this->get_ez_ext(0,0,0))*qm2 );

    //NOTE no cinv multiplied yet; see later
    bx0 = static_cast<real_long>( (bx[n] + this->get_bx_ext(0,0,0))*qm2 );
    by0 = static_cast<real_long>( (by[n] + this->get_by_ext(0,0,0))*qm2 );
    bz0 = static_cast<real_long>( (bz[n] + this->get_bz_ext(0,0,0))*qm2 );

    //-------------------------------------------------- 
    // first half electric acceleration
    u0 = c*vel0n + ex0;
    v0 = c*vel1n + ey0;
    w0 = c*vel2n + ez0;

    //-------------------------------------------------- 
    // intermediate gamma
    g2 = ( 1.0 + u0*u0 + v0*v0 + w0*w0 );
    b2 = bx0*bx0 + by0*by0 + bz0*bz0;

    // FIXME alternatively if cinv is in B_i terms
    //g2= (c*c + u0*u0 + v0*v0 + w0*w0)/(c*c);
    //b2 = (bx0*bx0 + by0*by0 + bz0*bz0)/c2;

    ginv = 1./sqrt( 0.5*(g2-b2 + sqrt( (g2-b2)*(g2-b2) + 4.0*(b2 + (bx0*u0 + by0*v0 + bz0*w0)*(bx0*u0 + by0*v0 + bz0*w0)))));

    //std::cout << "g_ vs gnew " << sqrt(g2) << " " << 1/ginv << "\n";

    //-------------------------------------------------- 
    // first half magnetic rotation
    bx0 *= ginv;
    by0 *= ginv;
    bz0 *= ginv;

    //std::cout << "bx prior rot" << bx0 << " " << by0 << " " << bz0 << "\n";

    f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
    u1 = (u0 + v0*bz0 - w0*by0)*f;
    v1 = (v0 + w0*bx0 - u0*bz0)*f;
    w1 = (w0 + u0*by0 - v0*bx0)*f;

    //-------------------------------------------------- 
    // second half of magnetic rotation & electric acceleration
    u0 = u0 + v1*bz0 - w1*by0 + ex0;
    v0 = v0 + w1*bx0 - u1*bz0 + ey0;
    w0 = w0 + u1*by0 - v1*bx0 + ez0;


    //-------------------------------------------------- 
    // normalized 4-velocity advance
    vel[0][n] = static_cast<real_prtcl>( u0*cinv );
    vel[1][n] = static_cast<real_prtcl>( v0*cinv );
    vel[2][n] = static_cast<real_prtcl>( w0*cinv );

    // position advance; 
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n]*ginv*c;
  }

}



//--------------------------------------------------
// explicit template instantiation

template class pic::HigueraCaryPusher<1,3>; // 1D3V
template class pic::HigueraCaryPusher<2,3>; // 2D3V
template class pic::HigueraCaryPusher<3,3>; // 3D3V

