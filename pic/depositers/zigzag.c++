#include "zigzag.h"

#include <algorithm>
#include <cassert>
#include <cmath>

using std::min;
using std::max;


// TODO: optimize for cases when we know D; now specialized for D=2
template<size_t D, size_t V>
void pic::ZigZag<D,V>::solve( pic::Tile<D>& tile )
{

  auto& yee = tile.get_yee();
  auto mins = tile.mins;

  for(auto&& container : tile.containers) {

    // initialize pointers to particle arrays
    int nparts = container.size();
      
    real_prtcl* loc[3];
    for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    real_prtcl* vel[3];
    for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    real_long invgam;
    real_long c = tile.cfl;
    real_long q = container.q;
    real_long x0, y0, z0, x1, x2, y1, y2, z1, z2;

    int i1,i2,j1,j2,k1,k2;

    real_long xr, yr, zr;
    real_long Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
    real_long Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;

    // loop and check particles
    int n1 = 0;
    int n2 = nparts;

    real_long loc0n, loc1n, loc2n, vel0n, vel1n, vel2n;


    // TODO: think SIMD (not possible due to ijk writing to yee)
    for(int n=n1; n<n2; n++) {

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

      // Fall into this debug section if something is very wrong; TODO: remove
      // NOTE: algorithm below has +1 on some indices so max ind is N + H-1 (hence 2 here)
      bool debug_flag = false;
      if (!( (i1 >= -3 && i1 < static_cast<int>(tile.mesh_lengths[0]+2)) )) debug_flag=true;
      if (!( (j1 >= -3 && j1 < static_cast<int>(tile.mesh_lengths[1]+2)) )) debug_flag=true;
      if (!( (k1 >= -3 && k1 < static_cast<int>(tile.mesh_lengths[2]+2)) )) debug_flag=true;

      if (!( (i2 >= -3 && i2 < static_cast<int>(tile.mesh_lengths[0]+2)) )) debug_flag=true;
      if (!( (j2 >= -3 && j2 < static_cast<int>(tile.mesh_lengths[1]+2)) )) debug_flag=true;
      if (!( (k2 >= -3 && k2 < static_cast<int>(tile.mesh_lengths[2]+2)) )) debug_flag=true;

      if (debug_flag) {
        std::cout << "--------------------------------------------------\n";
        std::cout << "n=" << n;
        std::cout << " i1: " << i1;
        std::cout << " j1: " << j1;
        std::cout << " k1: " << k1;
        std::cout << " ||| ";
        std::cout << " i2: " << i2;
        std::cout << " j2: " << j2;
        std::cout << " k2: " << k2;
        std::cout << "\n";

        std::cout << " x1sp: " << x1;
        std::cout << " y1sp: " << y1;
        std::cout << " z1sp: " << z1;
        std::cout << " ||| ";
        std::cout << " x2sp: " << x2;
        std::cout << " y2sp: " << y2;
        std::cout << " z2sp: " << z2;
        std::cout << " minxyz: " << mins[0] << " " << mins[1] << " " << mins[2];
        std::cout << " tilelen:" << tile.mesh_lengths[0] << " " << tile.mesh_lengths[1] << " " << tile.mesh_lengths[2];
        std::cout << "\n";

        std::cout << " vx: " <<  vel[0][n];
        std::cout << " vy: " <<  vel[1][n];
        std::cout << " vz: " <<  vel[2][n];
        std::cout << " gam: "<<  1.0/invgam;
        std::cout << "\n";

        std::cout << " xr: " <<  xr;
        std::cout << " yr: " <<  yr;
        std::cout << " zr: " <<  zr;
        std::cout << "\n";

        std::cout << " Fx1: " <<  Fx1;
        std::cout << " Fy1: " <<  Fy1;
        std::cout << " Fz1: " <<  Fz1;
        std::cout << " Wx1: " <<  Wx1;
        std::cout << " Wy1: " <<  Wy1;
        std::cout << " Wz1: " <<  Wz1;
        std::cout << "\n";

        std::cout << " Fx2: " <<  Fx2;
        std::cout << " Fy2: " <<  Fy2;
        std::cout << " Fz2: " <<  Fz2;
        std::cout << " Wx2: " <<  Wx2;
        std::cout << " Wy2: " <<  Wy2;
        std::cout << " Wz2: " <<  Wz2;
        std::cout << "\n";

        std::cout << std::flush;

        // always fail if we end here
        assert(true);
        continue;
      }

      //--------------------------------------------------
      // index checking
      // another check; TODO: remove
      //if(D >= 1 ) assert(i1 >= -3 && i1 < int(tile.mesh_lengths[0]+3)) ;
      //if(D >= 2 ) assert(j1 >= -3 && j1 < int(tile.mesh_lengths[1]+3)) ;
      //if(D >= 3 ) assert(k1 >= -3 && k1 < int(tile.mesh_lengths[2]+3)) ;
      //if (D >= 1) assert(i2 >= -3 && i2 < int(tile.mesh_lengths[0]+3));
      //if (D >= 2) assert(j2 >= -3 && j2 < int(tile.mesh_lengths[1]+3));
      //if (D >= 3) assert(k2 >= -3 && k2 < int(tile.mesh_lengths[2]+3));
        
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
      if(D>=1) yee.jz(i1+1, j1  , k1  ) += Fz1 * Wx1     * (1-Wy1);
      if(D>=2) yee.jz(i1  , j1+1, k1  ) += Fz1 * (1-Wx1) * Wy1;
      if(D>=2) yee.jz(i1+1, j1+1, k1  ) += Fz1 * Wx1     * Wy1;

      if(D>=1) yee.jz(i2  , j2  , k2  ) += Fz2 * (1-Wx2) * (1-Wy2);
      if(D>=1) yee.jz(i2+1, j2  , k2  ) += Fz2 * Wx2     * (1-Wy2);
      if(D>=2) yee.jz(i2  , j2+1, k2  ) += Fz2 * (1-Wx2) * Wy2;
      if(D>=2) yee.jz(i2+1, j2+1, k2  ) += Fz2 * Wx2     * Wy2;

    }

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag<1,3>; // 1D3V
template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag<3,3>; // 3D3V

