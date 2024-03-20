#include "boris_grav.h"

#include <cmath> 
#include "../../tools/signum.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using toolbox::sign;

template<size_t D, size_t V>
void pic::BorisPusherGrav<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const double c  = tile.cfl;
  const double qm = sign(con.q)/con.m; // q_s/m_s (sign only because emf are in units of q)
  const double m  = con.m; // mass
  

  // loop over particles
  UniIter::iterate([=] DEVCALLABLE (size_t n, pic::ParticleContainer<D>& con){

    double loc0n = con.loc(0,n);
    //double loc1n = con.loc(1,n);
    //double loc2n = con.loc(2,n);

    double vel0n = con.vel(0,n);
    double vel1n = con.vel(1,n);
    double vel2n = con.vel(2,n);

    // read particle-specific emf
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

    // addition of radiation pressure (gamma at half time step)
    // u at t + dt/2
    double uxt  = (u0/c + vel0n)*0.5;
    double uyt  = (v0/c + vel1n)*0.5;
    double uzt  = (w0/c + vel2n)*0.5;
    //ut   = sqrt(uxt*uxt + uyt*uyt + uzt*uzt);
    double gamt = sqrt(1.0 + uxt*uxt + uyt*uyt + uzt*uzt);
      
    //gravx = c*g0*gamt*(cenx - loc0n);
    double gravx = c*g0*m*gamt*sign(cenx - loc0n);
    double gravy = 0;
    double gravz = 0;

    //--------------------------------------------------
    // normalized 4-velocity advance
    con.vel(0,n) = u0/c + gravx;
    con.vel(1,n) = v0/c + gravy;
    con.vel(2,n) = w0/c + gravz;

    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) con.loc(i,n) += con.vel(i,n)*ginv*c;

  }, con.size(), con);

  UniIter::sync();


#ifdef GPU
  nvtxRangePop();
#endif
}


//--------------------------------------------------
// explicit template instantiation

template class pic::BorisPusherGrav<1,3>; // 1D3V
template class pic::BorisPusherGrav<2,3>; // 2D3V
template class pic::BorisPusherGrav<3,3>; // 3D3V
