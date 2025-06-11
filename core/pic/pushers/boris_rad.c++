#include <cmath> 

#include "core/pic/pushers/boris_rad.h"
#include "tools/signum.h"
#include "external/iter/iter.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::BorisPusherRad<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{

  const double c  = tile.cfl;
  const double qm = sign(con.q)/con.m; // q_s/m_s (sign only because emf are in units of q)


  // loop over particles
  UniIter::iterate([=]  (size_t n, pic::ParticleContainer<D>& con){
    
    // tmp variables
    double uxt, gamx, pressx, pressy, pressz;

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


    // rad pressure components; apply pressure only to particles behind the rad beam front
    if(loc0n <= beam_locx){

        // addition of radiation pressure (gamma at half time step)
        // u at t + dt/2
        uxt  = (u0/c + vel0n)*0.5;
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
    con.vel(0,n) = u0/c + pressx;
    con.vel(1,n) = v0/c + pressy;
    con.vel(2,n) = w0/c + pressz;


    // position advance
    // NOTE: no mixed-precision calc here. Can be problematic.
    ginv = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    for(size_t i=0; i<D; i++) con.loc(i,n) += con.vel(i,n)*ginv*c;

  }, con.size(), con);

}


//--------------------------------------------------
// explicit template instantiation

template class pic::BorisPusherRad<1,3>; // 1D3V
template class pic::BorisPusherRad<2,3>; // 2D3V
template class pic::BorisPusherRad<3,3>; // 3D3V
