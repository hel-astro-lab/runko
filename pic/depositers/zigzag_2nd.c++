#include "zigzag_2nd.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;


inline auto W2nd(float_m x, float_m xr, int i
        ) -> std::tuple<float_m,float_m,float_m>
{
    float_m xm = 0.5f*(x + xr) - i; // TODO: has +1 in it?

    float_m W1 = 0.125f*( 2.0f*( (xm+1.0f)-3.0f)*( (xm+1.0f)-3.0f) );
    float_m W2 = 0.75f- xm*xm;
    float_m W3 = 0.125f*( 2.0f*(-(xm-1.0f)-3.0f)*(-(xm-1.0f)-3.0f) );

    return { W1, W2, W3 };
}


template<size_t D, size_t V>
void pic::ZigZag_2nd<D,V>::solve( pic::Tile<D>& tile )
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


  for(auto&& con : tile.containers) {

    const double c = tile.cfl;    // speed of light
    const float_m q = con.q; // charge

    UniIter::iterate([=] DEVCALLABLE (
                size_t n, 
                fields::YeeLattice &yee,
                pic::ParticleContainer<D>& con
                ){

      //--------------------------------------------------
      double vel0n = con.vel(0,n);
      double vel1n = con.vel(1,n);
      double vel2n = con.vel(2,n);
      double invgam = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

      //--------------------------------------------------
      // new (normalized) location, x_{n+1}
        
      //float_m x1, y1, z1;
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

      auto [Wx1a, Wx2a, Wx3a] = W2nd(x1, xr, i1);
      auto [Wx1b, Wx2b, Wx3b] = W2nd(x2, xr, i2);

      auto [Wy1a, Wy2a, Wy3a] = W2nd(y1, yr, j1);
      auto [Wy1b, Wy2b, Wy3b] = W2nd(y2, yr, j2);

      auto [Wz1a, Wz2a, Wz3a] = W2nd(z1, zr, k1);
      auto [Wz1b, Wz2b, Wz3b] = W2nd(z2, zr, k2);

      //--------------------------------------------------
      //reduce dimension if needed
      if(D <= 1) {
        Wy1a, Wy2a, Wy3a = 1.0;
        Wy1b, Wy2b, Wy3b = 1.0;
      }

      if(D <= 2) {
        Wz1a, Wz2a, Wz3a = 1.0;
        Wz1b, Wz2b, Wz3b = 1.0;
      }

      // TODO need dimensionality switches to these too to optimize them
        
      const float_m Wxx1[3] = {Wx1a, Wx2a, Wx3a};
      const float_m Wxx2[3] = {Wx1a, Wx2a, Wx3a};

      const float_m Wyy1[3] = {Wy1a, Wy2a, Wy3a};
      const float_m Wyy2[3] = {Wy1a, Wy2a, Wy3a};

      const float_m Wzz1[3] = {Wz1a, Wz2a, Wz3a};
      const float_m Wzz2[3] = {Wz1a, Wz2a, Wz3a};

      //jx
      for(int zi=-1; zi <=1; zi++)
      for(int yi=-1; yi <=1; yi++){
        atomic_add( yee.jx(i1, j1+yi, k1+zi), Fx1*Wyy1[yi+1]*Wzz1[zi+1] );
        atomic_add( yee.jx(i2, j2+yi, k2+zi), Fx2*Wyy2[yi+1]*Wzz2[zi+1] );
      }

      //jy
      for(int zi=-1; zi <=1; zi++)
      for(int xi=-1; xi <=1; xi++){
        atomic_add( yee.jy(i1+xi, j1, k1+zi), Fy1*Wxx1[xi+1]*Wzz1[zi+1] );
        atomic_add( yee.jy(i2+xi, j2, k2+zi), Fy2*Wxx2[xi+1]*Wzz2[zi+1] );
      }

      //jz
      for(int yi=-1; yi <=1; yi++)
      for(int xi=-1; xi <=1; xi++){
        atomic_add( yee.jz(i1+xi, j1+yi, k1), Fz1*Wxx1[xi+1]*Wyy1[yi+1] );
        atomic_add( yee.jz(i2+xi, j2+yi, k2), Fz2*Wxx2[xi+1]*Wyy2[yi+1] );
      }


      //jx
      //atomic_add( yee.jx(i1  ,j1-1, k1-1), Fx1*Wy1a*Wz1a);
      //atomic_add( yee.jx(i1  ,j1  , k1-1), Fx1*Wy2a*Wz1a);
      //atomic_add( yee.jx(i1  ,j1+1, k1-1), Fx1*Wy3a*Wz1a);
      //atomic_add( yee.jx(i1  ,j1-1, k1  ), Fx1*Wy1a*Wz2a);
      //atomic_add( yee.jx(i1  ,j1  , k1  ), Fx1*Wy2a*Wz2a);
      //atomic_add( yee.jx(i1  ,j1+1, k1  ), Fx1*Wy3a*Wz2a);
      //atomic_add( yee.jx(i1  ,j1-1, k1+1), Fx1*Wy1a*Wz3a);
      //atomic_add( yee.jx(i1  ,j1  , k1+1), Fx1*Wy2a*Wz3a);
      //atomic_add( yee.jx(i1  ,j1+1, k1+1), Fx1*Wy3a*Wz3a);

      //atomic_add( yee.jx(i2  ,j2-1, k2-1), Fx2*Wy1b*Wz1b);
      //atomic_add( yee.jx(i2  ,j2  , k2-1), Fx2*Wy2b*Wz1b);
      //atomic_add( yee.jx(i2  ,j2+1, k2-1), Fx2*Wy3b*Wz1b);
      //atomic_add( yee.jx(i2  ,j2-1, k2  ), Fx2*Wy1b*Wz2b);
      //atomic_add( yee.jx(i2  ,j2  , k2  ), Fx2*Wy2b*Wz2b);
      //atomic_add( yee.jx(i2  ,j2+1, k2  ), Fx2*Wy3b*Wz2b);
      //atomic_add( yee.jx(i2  ,j2-1, k2+1), Fx2*Wy1b*Wz3b);
      //atomic_add( yee.jx(i2  ,j2  , k2+1), Fx2*Wy2b*Wz3b);
      //atomic_add( yee.jx(i2  ,j2+1, k2+1), Fx2*Wy3b*Wz3b);


      ////jy
      //atomic_add( yee.jy(i1-1,j1  , k1-1), Fy1*Wx1a*Wz1a);
      //atomic_add( yee.jy(i1  ,j1  , k1-1), Fy1*Wx2a*Wz1a);
      //atomic_add( yee.jy(i1+1,j1  , k1-1), Fy1*Wx3a*Wz1a);
      //atomic_add( yee.jy(i1-1,j1  , k1  ), Fy1*Wx1a*Wz2a);
      //atomic_add( yee.jy(i1  ,j1  , k1  ), Fy1*Wx2a*Wz2a);
      //atomic_add( yee.jy(i1+1,j1  , k1  ), Fy1*Wx3a*Wz2a);
      //atomic_add( yee.jy(i1-1,j1  , k1+1), Fy1*Wx1a*Wz3a);
      //atomic_add( yee.jy(i1  ,j1  , k1+1), Fy1*Wx2a*Wz3a);
      //atomic_add( yee.jy(i1+1,j1  , k1+1), Fy1*Wx3a*Wz3a);

      //atomic_add( yee.jy(i2-1,j2  , k2-1), Fy2*Wx1b*Wz1b);
      //atomic_add( yee.jy(i2  ,j2  , k2-1), Fy2*Wx2b*Wz1b);
      //atomic_add( yee.jy(i2+1,j2  , k2-1), Fy2*Wx3b*Wz1b);
      //atomic_add( yee.jy(i2-1,j2  , k2  ), Fy2*Wx1b*Wz2b);
      //atomic_add( yee.jy(i2  ,j2  , k2  ), Fy2*Wx2b*Wz2b);
      //atomic_add( yee.jy(i2+1,j2  , k2  ), Fy2*Wx3b*Wz2b);
      //atomic_add( yee.jy(i2-1,j2  , k2+1), Fy2*Wx1b*Wz3b);
      //atomic_add( yee.jy(i2  ,j2  , k2+1), Fy2*Wx2b*Wz3b);
      //atomic_add( yee.jy(i2+1,j2  , k2+1), Fy2*Wx3b*Wz3b);


      ////jz
      //atomic_add( yee.jy(i1-1,j1-1, k1  ), Fz1*Wx1a*Wy1a);
      //atomic_add( yee.jy(i1  ,j1-1, k1  ), Fz1*Wx2a*Wy1a);
      //atomic_add( yee.jy(i1+1,j1-1, k1  ), Fz1*Wx3a*Wy1a);
      //atomic_add( yee.jy(i1-1,j1  , k1  ), Fz1*Wx1a*Wy2a);
      //atomic_add( yee.jy(i1  ,j1  , k1  ), Fz1*Wx2a*Wy2a);
      //atomic_add( yee.jy(i1+1,j1  , k1  ), Fz1*Wx3a*Wy2a);
      //atomic_add( yee.jy(i1-1,j1+1, k1  ), Fz1*Wx1a*Wy3a);
      //atomic_add( yee.jy(i1  ,j1+1, k1  ), Fz1*Wx2a*Wy3a);
      //atomic_add( yee.jy(i1+1,j1+1, k1  ), Fz1*Wx3a*Wy3a);

      //atomic_add( yee.jy(i2-1,j2-1, k2  ), Fz2*Wx1b*Wy1b);
      //atomic_add( yee.jy(i2  ,j2-1, k2  ), Fz2*Wx2b*Wy1b);
      //atomic_add( yee.jy(i2+1,j2-1, k2  ), Fz2*Wx3b*Wy1b);
      //atomic_add( yee.jy(i2-1,j2  , k2  ), Fz2*Wx1b*Wy2b);
      //atomic_add( yee.jy(i2  ,j2  , k2  ), Fz2*Wx2b*Wy2b);
      //atomic_add( yee.jy(i2+1,j2  , k2  ), Fz2*Wx3b*Wy2b);
      //atomic_add( yee.jy(i2-1,j2+1, k2  ), Fz2*Wx1b*Wy3b);
      //atomic_add( yee.jy(i2  ,j2+1, k2  ), Fz2*Wx2b*Wy3b);
      //atomic_add( yee.jy(i2+1,j2+1, k2  ), Fz2*Wx3b*Wy3b);

    }, con.size(), yee, con);

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
//template class pic::ZigZag<1,3>; // 1D3V
//template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag_2nd<3,3>; // 3D3V

