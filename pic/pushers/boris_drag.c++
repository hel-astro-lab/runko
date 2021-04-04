#include "boris_drag.h"

#include <cmath> 
#include "../../tools/signum.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using toolbox::sign;


template<size_t D, size_t V>
double pic::BorisPusherDrag<D,V>::kn(double x)
{
  if (temp == 0.0) return 1.0; // Thomson scattering limit

  double sig;
  sig = (1.0 - 4.0/x - 8.0/x/x)*log(1.+x) + 0.5 + 8.0/x - 1.0/(2.0*(1. + x)*(1. + x));
  return (3.0/4.0)*sig/x;
}



template<size_t D, size_t V>
void pic::BorisPusherDrag<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{

  // maximum drag force experienced by particle
  const double dragthr = 0.1; 

  const double c  = tile.cfl;
  const double qm = sign(con.q)/con.m; // q_s/m_s (sign only because fields are in units of q)


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
    // Boris algorithm

    // first half electric acceleration
    double u0 = vel0n*c + ex0;
    double v0 = vel1n*c + ey0;
    double w0 = vel2n*c + ez0;

    // first half magnetic rotation
    double ginv = c/sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    bx0 *= ginv;
    by0 *= ginv;
    bz0 *= ginv;

    double f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
    double u1 = (u0 + v0*bz0 - w0*by0)*f;
    double v1 = (v0 + w0*bx0 - u0*bz0)*f;
    double w1 = (w0 + u0*by0 - v0*bx0)*f;

    // second half of magnetic rotation & electric acceleration
    u0 = u0 + v1*bz0 - w1*by0 + ex0;
    v0 = v0 + w1*bx0 - u1*bz0 + ey0;
    w0 = w0 + u1*by0 - v1*bx0 + ez0;


    //--------------------------------------------------
    // addition of drag (gamma at half time step)

    // u at t + dt/2
    double uxt  = (u0/c + vel0n)*0.5;
    double uyt  = (v0/c + vel1n)*0.5;
    double uzt  = (w0/c + vel2n)*0.5;
    double ut   = sqrt(uxt*uxt + uyt*uyt + uzt*uzt);
    double gamt = sqrt(1.0 + ut*ut);

    // subtract drag with Klein-Nishina reduction
    // A g^2 beta = A g^2 u/g = A g u
    double kncorr = kn(3.0*gamt*temp);

    // drag components
    double dragx = c*drag*kncorr*gamt*gamt*(uxt/gamt);
    double dragy = c*drag*kncorr*gamt*gamt*(uyt/gamt);
    double dragz = c*drag*kncorr*gamt*gamt*(uzt/gamt);

    // limit drag to maximum of dragthr of velocity
    double dragv = sqrt(dragx*dragx + dragy*dragy + dragz*dragz)/ut;
    double thr = 1.0;
    if (dragv > dragthr) thr = dragthr/dragv;

    //--------------------------------------------------
    // normalized 4-velocity advance with drag
    con.vel(0,n) = u0/c - thr*dragx;
    con.vel(1,n) = v0/c - thr*dragy;
    con.vel(2,n) = w0/c - thr*dragz;

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) con.loc(i,n) += con.vel(i,n)*ginv*c*freezing_factor;

  }, con.size(), con);

  UniIter::sync();


#ifdef GPU
  nvtxRangePop();
#endif
}



//--------------------------------------------------
// explicit template instantiation

template class pic::BorisPusherDrag<1,3>; // 1D3V
template class pic::BorisPusherDrag<2,3>; // 2D3V
template class pic::BorisPusherDrag<3,3>; // 3D3V
