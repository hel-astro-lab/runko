#include <cmath> 

#include "core/pic/pushers/boris.h"
#include "tools/signum.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


using toolbox::sign;

template<size_t D, size_t V>
void pic::BorisPusher<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const float c  = tile.cfl;
  const float qm = sign(con.q)/con.m; // q_s/m_s (sign only because fields are in units of q)

  // loop over particles
  UniIter::iterate([=] DEVCALLABLE (size_t n, pic::ParticleContainer<D>& con){
    float vel0n = con.vel(0,n)*c;
    float vel1n = con.vel(1,n)*c;
    float vel2n = con.vel(2,n)*c;

    // read particle-specific fields
    float ex0 = ( con.ex(n) + this->get_ex_ext(0,0,0) )*0.5*qm;
    float ey0 = ( con.ey(n) + this->get_ey_ext(0,0,0) )*0.5*qm;
    float ez0 = ( con.ez(n) + this->get_ez_ext(0,0,0) )*0.5*qm;

    float bx0 = ( con.bx(n) + this->get_bx_ext(0,0,0) )*0.5*qm/c;
    float by0 = ( con.by(n) + this->get_by_ext(0,0,0) )*0.5*qm/c;
    float bz0 = ( con.bz(n) + this->get_bz_ext(0,0,0) )*0.5*qm/c;

    //--------------------------------------------------
    // Boris algorithm

    // first half electric acceleration
    float u0 = vel0n + ex0;
    float v0 = vel1n + ey0;
    float w0 = vel2n + ez0;

    // first half magnetic rotation
    float ginv = c/sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
    bx0 *= ginv;
    by0 *= ginv;
    bz0 *= ginv;

    float f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
    float u1 = (u0 + v0*bz0 - w0*by0)*f;
    float v1 = (v0 + w0*bx0 - u0*bz0)*f;
    float w1 = (w0 + u0*by0 - v0*bx0)*f;

    // second half of magnetic rotation & electric acceleration
    u0 = u0 + v1*bz0 - w1*by0 + ex0;
    v0 = v0 + w1*bx0 - u1*bz0 + ey0;
    w0 = w0 + u1*by0 - v1*bx0 + ez0;


    //--------------------------------------------------
    // normalized 4-velocity advance
    con.vel(0,n) = u0/c;
    con.vel(1,n) = v0/c;
    con.vel(2,n) = w0/c;

    // position advance; 
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

template class pic::BorisPusher<1,3>; // 1D3V
template class pic::BorisPusher<2,3>; // 2D3V
template class pic::BorisPusher<3,3>; // 3D3V

