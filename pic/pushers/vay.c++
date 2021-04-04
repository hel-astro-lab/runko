#include "vay.h"

#include <cmath> 
#include "../../tools/signum.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using toolbox::sign;

template<size_t D, size_t V>
void pic::VayPusher<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{


#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const double c   = tile.cfl;
  const double cinv= 1.0/c;
  const double qm  = sign(con.q)/con.m; // q_s/m_s (sign only because fields are in units of q)
  

  // loop over particles
  UniIter::iterate([=] DEVCALLABLE (size_t n, pic::ParticleContainer<D>& con){

    double vel0n = con.vel(0,n);
    double vel1n = con.vel(1,n);
    double vel2n = con.vel(2,n);

    // read particle-specific fields
    double ex0 = ( con.ex(n) + this->get_ex_ext(0,0,0) )*0.5*qm;
    double ey0 = ( con.ey(n) + this->get_ey_ext(0,0,0) )*0.5*qm;
    double ez0 = ( con.ez(n) + this->get_ez_ext(0,0,0) )*0.5*qm;

    double bx0 = ( con.bx(n) + this->get_bx_ext(0,0,0) )*0.5*qm/c;
    double by0 = ( con.by(n) + this->get_by_ext(0,0,0) )*0.5*qm/c;
    double bz0 = ( con.bz(n) + this->get_bz_ext(0,0,0) )*0.5*qm/c;

    //--------------------------------------------------
    // Vay algorithm
      
    // gamma^-1
    double g = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);
    double vx0 = c*vel0n*g;
    double vy0 = c*vel1n*g;
    double vz0 = c*vel2n*g;

    // u' (cinv is already multiplied into B)
    double u1 = c*vel0n + 2.0*ex0 + vy0*bz0 - vz0*by0;
    double v1 = c*vel1n + 2.0*ey0 + vz0*bx0 - vx0*bz0;
    double w1 = c*vel2n + 2.0*ez0 + vx0*by0 - vy0*bx0;
    
    // gamma(u')
    double ustar = cinv*(u1*bx0+v1*by0+w1*bz0);
    double sig = cinv*cinv*( c*c + u1*u1+ v1*v1+ w1*w1) - (bx0*bx0+ by0*by0+ bz0*bz0);
    g = 1.0/sqrt( 0.5*(sig + sqrt(sig*sig + 4.0*(bx0*bx0 + by0*by0 + bz0*bz0 + ustar*ustar))));

    double tx = bx0*g;
    double ty = by0*g;
    double tz = bz0*g;
    double f = 1.0/(1.0+ tx*tx+ ty*ty+ tz*tz);
	  
    // final step 
    double u0 = f*(u1 + (u1*tx + v1*ty + w1*tz)*tx + v1*tz - w1*ty);
    double v0 = f*(v1 + (u1*tx + v1*ty + w1*tz)*ty + w1*tx - u1*tz);
    double w0 = f*(w1 + (u1*tx + v1*ty + w1*tz)*tz + u1*ty - v1*tx);


    //--------------------------------------------------
    // normalized 4-velocity advance
    con.vel(0,n) = u0/c;
    con.vel(1,n) = v0/c;
    con.vel(2,n) = w0/c;

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    g = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) con.loc(i,n) += con.vel(i,n)*g*c;

  }, con.size(), con);

  UniIter::sync();


#ifdef GPU
  nvtxRangePop();
#endif
}



//--------------------------------------------------
// explicit template instantiation

template class pic::VayPusher<1,3>; // 1D3V
template class pic::VayPusher<2,3>; // 2D3V
template class pic::VayPusher<3,3>; // 3D3V

