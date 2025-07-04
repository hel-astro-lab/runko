#include <cmath> 

#include "core/pic/pushers/higuera_cary.h"
#include "tools/signum.h"
#include "external/iter/iter.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::HigueraCaryPusher<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{

  const double c  = tile.cfl;
  const double qm = sign(con.q)/con.m; // q_s/m_s (sign only because emf are in units of q)


  // loop over particles
  UniIter::iterate([=]  (size_t n, pic::ParticleContainer<D>& con){
    double vel0n = con.vel(0,n);
    double vel1n = con.vel(1,n);
    double vel2n = con.vel(2,n);

    // read particle-specific emf
    double ex0 = ( con.ex(n) + this->get_ex_ext(0,0,0) )*0.5*qm;
    double ey0 = ( con.ey(n) + this->get_ey_ext(0,0,0) )*0.5*qm;
    double ez0 = ( con.ez(n) + this->get_ez_ext(0,0,0) )*0.5*qm;

    double bx0 = ( con.bx(n) + this->get_bx_ext(0,0,0) )*0.5*qm;
    double by0 = ( con.by(n) + this->get_by_ext(0,0,0) )*0.5*qm;
    double bz0 = ( con.bz(n) + this->get_bz_ext(0,0,0) )*0.5*qm;

    //-------------------------------------------------- 
    // first half electric acceleration
    double u0 = c*vel0n + ex0;
    double v0 = c*vel1n + ey0;
    double w0 = c*vel2n + ez0;

    //-------------------------------------------------- 
    // intermediate gamma
    double g2 = (c*c + u0*u0 + v0*v0 + w0*w0)/(c*c);
    double b2 = bx0*bx0 + by0*by0 + bz0*bz0;
    double ginv = 1./sqrt( 0.5*(g2-b2 + sqrt( (g2-b2)*(g2-b2) + 4.0*(b2 + (bx0*u0 + by0*v0 + bz0*w0)*(bx0*u0 + by0*v0 + bz0*w0)))));

    //-------------------------------------------------- 
    // first half magnetic rotation; cinv is multiplied to B field only here
    bx0 *= ginv/c;
    by0 *= ginv/c;
    bz0 *= ginv/c;

    double f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
    double u1 = (u0 + v0*bz0 - w0*by0)*f;
    double v1 = (v0 + w0*bx0 - u0*bz0)*f;
    double w1 = (w0 + u0*by0 - v0*bx0)*f;

    //-------------------------------------------------- 
    // second half of magnetic rotation & electric acceleration
    u0 = u0 + v1*bz0 - w1*by0 + ex0;
    v0 = v0 + w1*bx0 - u1*bz0 + ey0;
    w0 = w0 + u1*by0 - v1*bx0 + ez0;

    //-------------------------------------------------- 
    // normalized 4-velocity advance
    con.vel(0,n) = u0/c;
    con.vel(1,n) = v0/c;
    con.vel(2,n) = w0/c;

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) con.loc(i,n) += con.vel(i,n)*ginv*c;
  }, con.size(), con);

}


//--------------------------------------------------
// explicit template instantiation

template class pic::HigueraCaryPusher<1,3>; // 1D3V
template class pic::HigueraCaryPusher<2,3>; // 2D3V
template class pic::HigueraCaryPusher<3,3>; // 3D3V

