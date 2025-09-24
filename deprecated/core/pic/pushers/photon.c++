#include <cmath> 

#include "core/pic/pushers/photon.h"
#include "tools/signum.h"
#include "external/iter/iter.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::PhotonPusher<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile)
{
  const double c  = tile.cfl;

  // loop over particles
  UniIter::iterate([=]  (size_t n, pic::ParticleContainer<D>& con){
    double vel0n = con.vel(0,n);
    double vel1n = con.vel(1,n);
    double vel2n = con.vel(2,n);

    // normalize
    double u = sqrt(vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

    // first half electric acceleration
    double u0 = vel0n/u;
    double v0 = vel1n/u;
    double w0 = vel2n/u;

    //std::cout << "pushing " 
    //<< n << " u:" 
    //<< u0 << " " << v0 << " " << w0 << " loc:"
    //<< con.loc(0,n) << " " << con.loc(1,n) << " " << con.loc(2,n) << "\n";

    // position advance
    if(D >= 1) con.loc(0,n) += u0*c;
    if(D >= 2) con.loc(1,n) += v0*c;
    if(D >= 3) con.loc(2,n) += w0*c;

  }, con.size(), con);

}



//--------------------------------------------------
// explicit template instantiation

template class pic::PhotonPusher<1,3>; // 1D3V
template class pic::PhotonPusher<2,3>; // 2D3V
template class pic::PhotonPusher<3,3>; // 3D3V

