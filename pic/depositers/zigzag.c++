#include "zigzag.h"

#include <algorithm>
#include <cassert>
#include <cmath>

using std::min;
using std::max;

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// vectorization is broken at the point where we add values back to the grid.
// This auxiliary function tries to hide that complexity.
template<typename T, typename S>
DEVCALLABLE inline void atomic_add(T& lhs, S rhs) 
{
#ifdef GPU
    atomicAdd(&lhs, static_cast<T>(rhs));
#else
    //NOTE: need to use #pragma omp atomic if vectorizing these
#pragma omp atomic update
    lhs += static_cast<S>(rhs);
#endif
}



template<size_t D, size_t V>
void pic::ZigZag<D,V>::solve( pic::Tile<D>& tile )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& yee = tile.get_yee();
  const auto mins = tile.mins;

  //clear arrays before new update
  yee.jx.clear();
  yee.jy.clear();
  yee.jz.clear();
  yee.rho.clear();


  for(auto&& con: tile.containers) {

    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge

//    // no vectorization here since we dont use the general iterator
//    for(size_t n=0; n<con.size(); n++) {

    UniIter::iterate([=] DEVCALLABLE (
                size_t n, 
                fields::YeeLattice &yee,
                pic::ParticleContainer<D>& con
                ){

      //--------------------------------------------------
      // NOTE: performing velocity calculations via doubles to retain accuracy
      double vel0n = con.vel(0,n);
      double vel1n = con.vel(1,n);
      double vel2n = con.vel(2,n);
      double invgam = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

      //--------------------------------------------------
      // new (normalized) location, x_{n+1}
      float_m x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      float_m y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      float_m z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      // previos location, x_n
      float_m x1 = x2 - vel0n*invgam*c;
      float_m y1 = y2 - vel1n*invgam*c;
      float_m z1 = z2 - vel2n*invgam*c; 

      //--------------------------------------------------
      int i1  = D >= 1 ? floor(x1) : 0;
      int i2  = D >= 1 ? floor(x2) : 0;
      int j1  = D >= 2 ? floor(y1) : 0;
      int j2  = D >= 2 ? floor(y2) : 0;
      int k1  = D >= 3 ? floor(z1) : 0;
      int k2  = D >= 3 ? floor(z2) : 0;

      // relay point; +1 is equal to +\Delta x
      float_m xr = min( float_m(min(i1,i2)+1), max( float_m(max(i1,i2)), float_m(0.5*(x1+x2)) ) );
      float_m yr = min( float_m(min(j1,j2)+1), max( float_m(max(j1,j2)), float_m(0.5*(y1+y2)) ) );
      float_m zr = min( float_m(min(k1,k2)+1), max( float_m(max(k1,k2)), float_m(0.5*(z1+z2)) ) );

      //--------------------------------------------------
      // +q since - sign is already included in the Ampere's equation
      //q = weight*qe;
      float_m Fx1 = +q*(xr - x1);
      float_m Fy1 = +q*(yr - y1);
      float_m Fz1 = +q*(zr - z1);
      
      float_m Fx2 = +q*(x2 - xr);
      float_m Fy2 = +q*(y2 - yr);
      float_m Fz2 = +q*(z2 - zr);


      float_m Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
      float_m Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
      float_m Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

      float_m Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
      float_m Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
      float_m Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;


      //--------------------------------------------------
      // jx
      if(D>=1) atomic_add( yee.jx(i1  , j1  , k1  ), Fx1*(1.0f-Wy1)*(1.0f-Wz1) );
      if(D>=2) atomic_add( yee.jx(i1  , j1+1, k1  ), Fx1*Wy1       *(1.0f-Wz1) );
      if(D>=3) atomic_add( yee.jx(i1  , j1  , k1+1), Fx1*(1.0f-Wy1)*Wz1        );
      if(D>=3) atomic_add( yee.jx(i1  , j1+1, k1+1), Fx1*Wy1       *Wz1        );

      if(D>=1) atomic_add( yee.jx(i2  , j2  , k2  ), Fx2*(1.0f-Wy2)*(1.0f-Wz2) );
      if(D>=2) atomic_add( yee.jx(i2  , j2+1, k2  ), Fx2*Wy2       *(1.0f-Wz2) );
      if(D>=3) atomic_add( yee.jx(i2  , j2  , k2+1), Fx2*(1.0f-Wy2)*Wz2        );
      if(D>=3) atomic_add( yee.jx(i2  , j2+1, k2+1), Fx2*Wy2       *Wz2        );

      // jy
      if(D>=1) atomic_add( yee.jy(i1  , j1  , k1  ), Fy1*(1.0f-Wx1)*(1.0f-Wz1) );
      if(D>=2) atomic_add( yee.jy(i1+1, j1  , k1  ), Fy1*Wx1       *(1.0f-Wz1) );
      if(D>=3) atomic_add( yee.jy(i1  , j1  , k1+1), Fy1*(1.0f-Wx1)*Wz1        );
      if(D>=3) atomic_add( yee.jy(i1+1, j1  , k1+1), Fy1*Wx1       *Wz1        );

      if(D>=1) atomic_add( yee.jy(i2  , j2  , k2  ), Fy2*(1.0f-Wx2)*(1.0f-Wz2) );
      if(D>=2) atomic_add( yee.jy(i2+1, j2  , k2  ), Fy2*Wx2       *(1.0f-Wz2) );
      if(D>=3) atomic_add( yee.jy(i2  , j2  , k2+1), Fy2*(1.0f-Wx2)*Wz2        );
      if(D>=3) atomic_add( yee.jy(i2+1, j2  , k2+1), Fy2*Wx2       *Wz2        );

      // jz
      if(D>=1) atomic_add( yee.jz(i1  , j1  , k1  ), Fz1*(1.0f-Wx1)*(1.0f-Wy1) );
      if(D>=2) atomic_add( yee.jz(i1+1, j1  , k1  ), Fz1*Wx1       *(1.0f-Wy1) );
      if(D>=3) atomic_add( yee.jz(i1  , j1+1, k1  ), Fz1*(1.0f-Wx1)*Wy1        );
      if(D>=3) atomic_add( yee.jz(i1+1, j1+1, k1  ), Fz1*Wx1       *Wy1        );

      if(D>=1) atomic_add( yee.jz(i2  , j2  , k2  ), Fz2*(1.0f-Wx2)*(1.0f-Wy2) );
      if(D>=1) atomic_add( yee.jz(i2+1, j2  , k2  ), Fz2*Wx2       *(1.0f-Wy2) );
      if(D>=2) atomic_add( yee.jz(i2  , j2+1, k2  ), Fz2*(1.0f-Wx2)*Wy2        );
      if(D>=2) atomic_add( yee.jz(i2+1, j2+1, k2  ), Fz2*Wx2       *Wy2        );

    }, con.size(), yee, con);

    UniIter::sync();
  }//end of loop over species


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag<1,3>; // 1D3V
template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag<3,3>; // 3D3V
