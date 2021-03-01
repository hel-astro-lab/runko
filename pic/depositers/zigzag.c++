#include "zigzag.h"

#include <algorithm>
#include <cassert>
#include <cmath>

using std::min;
using std::max;

#include <nvtx3/nvToolsExt.h> 

#include "../../tools/iter/iter.h"


template<size_t D, size_t V>
void pic::ZigZag<D,V>::solve( pic::Tile<D>& tile )
{

  nvtxRangePush(__PRETTY_FUNCTION__);

  auto& yee = tile.get_yee();
  yee.jx.clear();
  yee.jy.clear();
  yee.jz.clear();
  yee.rho.clear();

  auto mins = tile.mins;

  for(auto&& container : tile.containers) {

    // initialize pointers to particle arrays
    int nparts = container.size();
      
    real_prtcl* loc[3];
    for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    real_prtcl* vel[3];
    for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );
    // loop and check particles
    int n1 = 0;
    int n2 = nparts;

    real_long c = tile.cfl;
    real_long q = container.q;

    #ifdef GPU
    UniIter::UniIterCU::iterate([=] __device__ (int n, fields::YeeLattice &yee){
      real_long invgam;
      real_long x0, y0, z0, x1, x2, y1, y2, z1, z2;

      int i1,i2,j1,j2,k1,k2;

      real_long xr, yr, zr;
      real_long Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
      real_long Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;


      real_long loc0n, loc1n, loc2n, vel0n, vel1n, vel2n;


      loc0n = static_cast<real_long>(loc[0][n]);
      loc1n = static_cast<real_long>(loc[1][n]);
      loc2n = static_cast<real_long>(loc[2][n]);

      vel0n = static_cast<real_long>(vel[0][n]);
      vel1n = static_cast<real_long>(vel[1][n]);
      vel2n = static_cast<real_long>(vel[2][n]);

      invgam = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

      x0 = loc0n - vel0n*invgam*c;
      y0 = loc1n - vel1n*invgam*c;
      z0 = loc2n - vel2n*invgam*c; 

      // normalized location w.r.t. tile; previous loc (x1) and current loc (x2)
      x1 = D >= 1 ? x0     - mins[0] : x0;
      x2 = D >= 1 ? loc0n  - mins[0] : loc0n;
      y1 = D >= 2 ? y0     - mins[1] : y0;
      y2 = D >= 2 ? loc1n  - mins[1] : loc1n;
      z1 = D >= 3 ? z0     - mins[2] : z0;
      z2 = D >= 3 ? loc2n  - mins[2] : loc2n;

  	  i1  = D >= 1 ? static_cast<int>(floor( x1 )) : 0;
  	  i2  = D >= 1 ? static_cast<int>(floor( x2 )) : 0;
  	  j1  = D >= 2 ? static_cast<int>(floor( y1 )) : 0;
  	  j2  = D >= 2 ? static_cast<int>(floor( y2 )) : 0;
  	  k1  = D >= 3 ? static_cast<int>(floor( z1 )) : 0;
  	  k2  = D >= 3 ? static_cast<int>(floor( z2 )) : 0;

      // relay point; +1 is equal to +\Delta x
      xr = min( real_long(min(i1,i2)+1), max( (real_long)max(i1,i2), (real_long)0.5*(x1+x2) ) );
      yr = min( real_long(min(j1,j2)+1), max( (real_long)max(j1,j2), (real_long)0.5*(y1+y2) ) );
      zr = min( real_long(min(k1,k2)+1), max( (real_long)max(k1,k2), (real_long)0.5*(z1+z2) ) );

      // +q since - sign is already included in the Ampere's equation
      //q = weight*qe;
      Fx1 = +q*(xr - x1);
      Fy1 = +q*(yr - y1);
      Fz1 = +q*(zr - z1);
      
      Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
      Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
      Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

      Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
      Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  	  Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;

      Fx2 = +q*(x2-xr);
      Fy2 = +q*(y2-yr);
      Fz2 = +q*(z2-zr);

      //--------------------------------------------------
      // index checking

      // jx
      atomicAdd(&yee.jx(i1,  j1,   k1)  , Fx1 * (1.0-Wy1) * (1.0-Wz1));
      atomicAdd(&yee.jx(i1,  j1+1, k1)  , Fx1 * Wy1       * (1.0-Wz1));
      atomicAdd(&yee.jx(i1,  j1,   k1+1), Fx1 * (1.0-Wy1) * Wz1);
      atomicAdd(&yee.jx(i1,  j1+1, k1+1), Fx1 * Wy1       * Wz1);

      atomicAdd(&yee.jx(i2,  j2,   k2)  , Fx2 * (1.0-Wy2) * (1.0-Wz2));
      atomicAdd(&yee.jx(i2,  j2+1, k2)  , Fx2 * Wy2       * (1.0-Wz2));
      atomicAdd(&yee.jx(i2,  j2,   k2+1), Fx2 * (1.0-Wy2) * Wz2);
      atomicAdd(&yee.jx(i2,  j2+1, k2+1), Fx2 * Wy2       * Wz2);

      // jy
      atomicAdd(&yee.jy(i1,  j1,   k1)  , Fy1 * (1.0-Wx1) * (1.0-Wz1));
      atomicAdd(&yee.jy(i1+1,j1,   k1)  , Fy1 * Wx1       * (1.0-Wz1));
      atomicAdd(&yee.jy(i1  ,j1,   k1+1), Fy1 * (1.0-Wx1) * Wz1);
      atomicAdd(&yee.jy(i1+1,j1,   k1+1), Fy1 * Wx1       * Wz1);
      
      atomicAdd(&yee.jy(i2,  j2,   k2)  , Fy2 * (1.0-Wx2) * (1.0-Wz2));
      atomicAdd(&yee.jy(i2+1,j2,   k2)  , Fy2 * Wx2       * (1.0-Wz2));
      atomicAdd(&yee.jy(i2,  j2,   k2+1), Fy2 * (1.0-Wx2) * Wz2);
      atomicAdd(&yee.jy(i2+1,j2,   k2+1), Fy2 * Wx2       * Wz2);
                            
      // jz
      atomicAdd(&yee.jz(i1,  j1,   k1)  , Fz1 * (1.0-Wx1) * (1.0-Wy1));
      atomicAdd(&yee.jz(i1+1,j1,   k1)  , Fz1 * Wx1       * (1.0-Wy1));
      atomicAdd(&yee.jz(i1,  j1+1, k1)  , Fz1 * (1.0-Wx1) * Wy1);
      atomicAdd(&yee.jz(i1+1,j1+1, k1)  , Fz1 * Wx1       * Wy1);

      atomicAdd(&yee.jz(i2,  j2,   k2)  , Fz2 * (1.0-Wx2) * (1.0-Wy2));
      atomicAdd(&yee.jz(i2+1,j2,   k2)  , Fz2 * Wx2       * (1.0-Wy2));
      atomicAdd(&yee.jz(i2,  j2+1, k2)  , Fz2 * (1.0-Wx2) * Wy2);
      atomicAdd(&yee.jz(i2+1,j2+1, k2)  , Fz2 * Wx2       * Wy2);

    }, nparts, yee);

    UniIter::sync();

    #else
    // TODO: think SIMD (not possible due to ijk writing to yee)
    //#pragma omp parallel for
    for(int n=n1; n<n2; n++) {
      real_long invgam;
      real_long x0, y0, z0, x1, x2, y1, y2, z1, z2;

      int i1,i2,j1,j2,k1,k2;

      real_long xr, yr, zr;
      real_long Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
      real_long Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;


      real_long loc0n, loc1n, loc2n, vel0n, vel1n, vel2n;


      loc0n = static_cast<real_long>(loc[0][n]);
      loc1n = static_cast<real_long>(loc[1][n]);
      loc2n = static_cast<real_long>(loc[2][n]);

      vel0n = static_cast<real_long>(vel[0][n]);
      vel1n = static_cast<real_long>(vel[1][n]);
      vel2n = static_cast<real_long>(vel[2][n]);

      invgam = 1.0/sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n);

      x0 = loc0n - vel0n*invgam*c;
      y0 = loc1n - vel1n*invgam*c;
      z0 = loc2n - vel2n*invgam*c; 

      // normalized location w.r.t. tile; previous loc (x1) and current loc (x2)
      // NOTE: these are floored because the prtcl orginates from inside the grid
      //       and that mapping is done via flooring.
      x1 = D >= 1 ? x0 - mins[0] : x0;
      y1 = D >= 2 ? y0 - mins[1] : y0;
      z1 = D >= 3 ? z0 - mins[2] : z0;
      i1 = D >= 1 ? static_cast<int>( floor(x1) ) : 0;
      j1 = D >= 2 ? static_cast<int>( floor(y1) ) : 0;
      k1 = D >= 3 ? static_cast<int>( floor(z1) ) : 0;

      x2 = D >= 1 ? loc0n  - mins[0] : loc0n;
      y2 = D >= 2 ? loc1n  - mins[1] : loc1n;
      z2 = D >= 3 ? loc2n  - mins[2] : loc2n;
      i2 = D >= 1 ? static_cast<int>( floor(x2) ) : 0;
      j2 = D >= 2 ? static_cast<int>( floor(y2) ) : 0;
      k2 = D >= 3 ? static_cast<int>( floor(z2) ) : 0;

      //// relay point; +1 is equal to +\Delta x
      xr = min( real_long(min(i1,i2)+1), max( real_long(max(i1,i2)), real_long(0.5*(x1+x2)) ) );
      yr = min( real_long(min(j1,j2)+1), max( real_long(max(j1,j2)), real_long(0.5*(y1+y2)) ) );
      zr = min( real_long(min(k1,k2)+1), max( real_long(max(k1,k2)), real_long(0.5*(z1+z2)) ) );

      // +q since - sign is already included in the Ampere's equation
      //q = weight*qe;
      Fx1 = +q*(xr - x1);
      Fy1 = +q*(yr - y1);
      Fz1 = +q*(zr - z1);
      
      Wx1 = D >= 1 ? 0.5*(x1 + xr) - i1 : 0.0;
      Wy1 = D >= 2 ? 0.5*(y1 + yr) - j1 : 0.0;
      Wz1 = D >= 3 ? 0.5*(z1 + zr) - k1 : 0.0;

      Wx2 = D >= 1 ? 0.5*(x2 + xr) - i2 : 0.0;
      Wy2 = D >= 2 ? 0.5*(y2 + yr) - j2 : 0.0;
  	  Wz2 = D >= 3 ? 0.5*(z2 + zr) - k2 : 0.0;

      Fx2 = +q*(x2-xr);
      Fy2 = +q*(y2-yr);
      Fz2 = +q*(z2-zr);

        
      // jx
      if(D>=1) yee.jx(i1  , j1  , k1  ) += Fx1 * (1-Wy1) * (1-Wz1);
      if(D>=2) yee.jx(i1  , j1+1, k1  ) += Fx1 * Wy1     * (1-Wz1);
      if(D>=3) yee.jx(i1  , j1  , k1+1) += Fx1 * (1-Wy1) * Wz1;
      if(D>=3) yee.jx(i1  , j1+1, k1+1) += Fx1 * Wy1     * Wz1;

      if(D>=1) yee.jx(i2  , j2  , k2  ) += Fx2 * (1-Wy2) * (1-Wz2);
      if(D>=2) yee.jx(i2  , j2+1, k2  ) += Fx2 * Wy2     * (1-Wz2);
      if(D>=3) yee.jx(i2  , j2  , k2+1) += Fx2 * (1-Wy2) * Wz2;
      if(D>=3) yee.jx(i2  , j2+1, k2+1) += Fx2 * Wy2     * Wz2;

      //// jy
      if(D>=1) yee.jy(i1  , j1  , k1  ) += Fy1 * (1-Wx1) * (1-Wz1);
      if(D>=2) yee.jy(i1+1, j1  , k1  ) += Fy1 * Wx1     * (1-Wz1);
      if(D>=3) yee.jy(i1  , j1  , k1+1) += Fy1 * (1-Wx1) * Wz1;
      if(D>=3) yee.jy(i1+1, j1  , k1+1) += Fy1 * Wx1     * Wz1;

      if(D>=1) yee.jy(i2  , j2  , k2  ) += Fy2 * (1-Wx2) * (1-Wz2);
      if(D>=2) yee.jy(i2+1, j2  , k2  ) += Fy2 * Wx2     * (1-Wz2);
      if(D>=3) yee.jy(i2  , j2  , k2+1) += Fy2 * (1-Wx2) * Wz2;
      if(D>=3) yee.jy(i2+1, j2  , k2+1) += Fy2 * Wx2     * Wz2;

      //// jz
      if(D>=1) yee.jz(i1  , j1  , k1  ) += Fz1 * (1-Wx1) * (1-Wy1);
      if(D>=2) yee.jz(i1+1, j1  , k1  ) += Fz1 * Wx1     * (1-Wy1);
      if(D>=3) yee.jz(i1  , j1+1, k1  ) += Fz1 * (1-Wx1) * Wy1;
      if(D>=3) yee.jz(i1+1, j1+1, k1  ) += Fz1 * Wx1     * Wy1;

      if(D>=1) yee.jz(i2  , j2  , k2  ) += Fz2 * (1-Wx2) * (1-Wy2);
      if(D>=1) yee.jz(i2+1, j2  , k2  ) += Fz2 * Wx2     * (1-Wy2);
      if(D>=1) yee.jz(i2  , j2+1, k2  ) += Fz2 * (1-Wx2) * Wy2;
      if(D>=1) yee.jz(i2+1, j2+1, k2  ) += Fz2 * Wx2     * Wy2;

      //#pragma omp atomic
      yee.jx(i1,  j1,   k1)   += Fx1 * (1.0-Wy1) * (1.0-Wz1);
      //#pragma omp atomic
      yee.jx(i1,  j1+1, k1)   += Fx1 * Wy1       * (1.0-Wz1);
      //#pragma omp atomic
      yee.jx(i1,  j1,   k1+1) += Fx1 * (1.0-Wy1) * Wz1;
      //#pragma omp atomic
      yee.jx(i1,  j1+1, k1+1) += Fx1 * Wy1       * Wz1;

      //#pragma omp atomic
      yee.jx(i2,  j2,   k2)   += Fx2 * (1.0-Wy2) * (1.0-Wz2);
      //#pragma omp atomic
      yee.jx(i2,  j2+1, k2)   += Fx2 * Wy2       * (1.0-Wz2);
      //#pragma omp atomic
      yee.jx(i2,  j2,   k2+1) += Fx2 * (1.0-Wy2) * Wz2;
      //#pragma omp atomic
      yee.jx(i2,  j2+1, k2+1) += Fx2 * Wy2       * Wz2;

      // jy
      //#pragma omp atomic
      yee.jy(i1,  j1,   k1)   += Fy1 * (1.0-Wx1) * (1.0-Wz1);
      //#pragma omp atomic
      yee.jy(i1+1,j1,   k1)   += Fy1 * Wx1       * (1.0-Wz1);
      //#pragma omp atomic
      yee.jy(i1  ,j1,   k1+1) += Fy1 * (1.0-Wx1) * Wz1;
      //#pragma omp atomic
      yee.jy(i1+1,j1,   k1+1) += Fy1 * Wx1       * Wz1;
      
      //#pragma omp atomic
      yee.jy(i2,  j2,   k2)   += Fy2 * (1.0-Wx2) * (1.0-Wz2);
      //#pragma omp atomic
      yee.jy(i2+1,j2,   k2)   += Fy2 * Wx2       * (1.0-Wz2);
      //#pragma omp atomic
      yee.jy(i2,  j2,   k2+1) += Fy2 * (1.0-Wx2) * Wz2;
      //#pragma omp atomic
      yee.jy(i2+1,j2,   k2+1) += Fy2 * Wx2       * Wz2;
                            
      // jz
      //#pragma omp atomic
      yee.jz(i1,  j1,   k1)   += Fz1 * (1.0-Wx1) * (1.0-Wy1);
      //#pragma omp atomic
      yee.jz(i1+1,j1,   k1)   += Fz1 * Wx1       * (1.0-Wy1);
      //#pragma omp atomic
      yee.jz(i1,  j1+1, k1)   += Fz1 * (1.0-Wx1) * Wy1;
      //#pragma omp atomic
      yee.jz(i1+1,j1+1, k1)   += Fz1 * Wx1       * Wy1;

      //#pragma omp atomic
      yee.jz(i2,  j2,   k2)   += Fz2 * (1.0-Wx2) * (1.0-Wy2);
      //#pragma omp atomic
      yee.jz(i2+1,j2,   k2)   += Fz2 * Wx2       * (1.0-Wy2);
      //#pragma omp atomic
      yee.jz(i2,  j2+1, k2)   += Fz2 * (1.0-Wx2) * Wy2;
      //#pragma omp atomic
      yee.jz(i2+1,j2+1, k2)   += Fz2 * Wx2       * Wy2;

    }
    #endif

  }//end of loop over species

  nvtxRangePop();
}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag<1,3>; // 1D3V
template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag<3,3>; // 3D3V

