#include "zigzag.h"

#include <algorithm>
#include <cassert>
#include <math.h>

using std::min;
using std::max;


// TODO: optimize for cases when we know D; now specialized for D=2
template<size_t D, size_t V>
void pic::ZigZag<D,V>::solve( pic::Tile<D>& tile )
{

  auto& yee = tile.get_yee();
  yee.jx.clear();
  yee.jy.clear();
  yee.jz.clear();


  auto mins = tile.mins;

  for(auto&& container : tile.containers) {

    // initialize pointers to particle arrays
    int nparts = container.size();
      
    double* loc[3];
    for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    double invgam;
    double c = tile.cfl;
    double q = container.q;
    // TODO: remove
    //std::cout << " q = " << q << " ispc: " << ispc << '\n';

    double x0, y0, z0, x1, x2, y1, y2, z1, z2;

    int i1,i2,j1,j2,k1,k2;

    double xr, yr, zr;

    double Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
    double Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;

    // loop and check particles
    int n1 = 0;
    int n2 = nparts;


    // TODO: think SIMD (not possible due to ijk writing to yee)
    for(int n=n1; n<n2; n++) {

      invgam = 1.0/sqrt(1.0 + 
          vel[0][n]*vel[0][n] + 
          vel[1][n]*vel[1][n] + 
          vel[2][n]*vel[2][n]);

      x0 = loc[0][n] - vel[0][n]*invgam*c;
      y0 = loc[1][n] - vel[1][n]*invgam*c;
      z0 = loc[2][n] - vel[2][n]*invgam*c; 


      // normalized location w.r.t. tile
      x1 = D >= 1 ? x0         - mins[0] : x0;
      x2 = D >= 1 ? loc[0][n]  - mins[0] : loc[0][n];
      y1 = D >= 2 ? y0         - mins[1] : y0;
      y2 = D >= 2 ? loc[1][n]  - mins[1] : loc[1][n];
      z1 = D >= 3 ? z0         - mins[2] : z0;
      z2 = D >= 3 ? loc[2][n]  - mins[2] : loc[2][n];

      // XXX: original algorithm has int(float(x)) != floor(x)
  	  i1  = D >= 1 ? static_cast<int>(floor( x1 )) : 0;
  	  i2  = D >= 1 ? static_cast<int>(floor( x2 )) : 0;
  	  j1  = D >= 2 ? static_cast<int>(floor( y1 )) : 0;
  	  j2  = D >= 2 ? static_cast<int>(floor( y2 )) : 0;
  	  k1  = D >= 3 ? static_cast<int>(floor( z1 )) : 0;
  	  k2  = D >= 3 ? static_cast<int>(floor( z2 )) : 0;
  	  //i1  = D >= 1 ? static_cast<int>(x1) : 0;
  	  //i2  = D >= 1 ? static_cast<int>(x2) : 0;
  	  //j1  = D >= 2 ? static_cast<int>(y1) : 0;
  	  //j2  = D >= 2 ? static_cast<int>(y2) : 0;
  	  //k1  = D >= 3 ? static_cast<int>(z1) : 0;
  	  //k2  = D >= 3 ? static_cast<int>(z2) : 0;

      // relay point; +1 is equal to +\Delta x
      xr = min( double(min(i1,i2)+1), max( double(max(i1,i2)), 0.5*(double(x1+x2)) ) );
      yr = min( double(min(j1,j2)+1), max( double(max(j1,j2)), 0.5*(double(y1+y2)) ) );
      zr = min( double(min(k1,k2)+1), max( double(max(k1,k2)), 0.5*(double(z1+z2)) ) );

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
      bool debug_flag = false;
      if (!( (i1 >= -3 && i1 < static_cast<int>(tile.mesh_lengths[0]+3)) )) debug_flag=true;
      if (!( (j1 >= -3 && j1 < static_cast<int>(tile.mesh_lengths[1]+3)) )) debug_flag=true;
      if (!( (k1 >= -3 && k1 < static_cast<int>(tile.mesh_lengths[2]+3)) )) debug_flag=true;

      if (!( (i2 >= -3 && i2 < static_cast<int>(tile.mesh_lengths[0]+3)) )) debug_flag=true;
      if (!( (j2 >= -3 && j2 < static_cast<int>(tile.mesh_lengths[1]+3)) )) debug_flag=true;
      if (!( (k2 >= -3 && k2 < static_cast<int>(tile.mesh_lengths[2]+3)) )) debug_flag=true;

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
        std::cout << "\n";

        std::cout << " vx: " <<  vel[0][n];
        std::cout << " vy: " <<  vel[1][n];
        std::cout << " vz: " <<  vel[2][n];
        std::cout << " gam: "<<  invgam;
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
        assert(false);
      }

      //--------------------------------------------------
      // index checking
      // another check; TODO: remove
      if(D >= 1 ) assert(i1   >= -3 && i1   < int(tile.mesh_lengths[0]+3)) ;
      if(D >= 2 ) assert(j1   >= -3 && j1   < int(tile.mesh_lengths[1]+3)) ;
      if(D >= 3 ) assert(k1   >= -3 && k1   < int(tile.mesh_lengths[2]+3)) ;

      if (D >= 1) assert(i2   >= -3 && i2   < int(tile.mesh_lengths[0]+3));
      if (D >= 2) assert(j2   >= -3 && j2   < int(tile.mesh_lengths[1]+3));
      if (D >= 3) assert(k2   >= -3 && k2   < int(tile.mesh_lengths[2]+3));


      // jx
      if (D >= 1) yee.jx(i1,  j1,   k1)   += Fx1 * (1.0-Wy1) * (1.0-Wz1);
      if (D >= 2) yee.jx(i1,  j1+1, k1)   += Fx1 * Wy1       * (1.0-Wz1);
      if (D >= 3) yee.jx(i1,  j1,   k1+1) += Fx1 * (1.0-Wy1) * Wz1;
      if (D >= 3) yee.jx(i1,  j1+1, k1+1) += Fx1 * Wy1       * Wz1;

      if (D >= 1) yee.jx(i2,  j2,   k2)   += Fx2 * (1.0-Wy2) * (1.0-Wz2);
      if (D >= 2) yee.jx(i2,  j2+1, k2)   += Fx2 * Wy2       * (1.0-Wz2);
      if (D >= 3) yee.jx(i2,  j2,   k2+1) += Fx2 * (1.0-Wy2) * Wz2;
      if (D >= 3) yee.jx(i2,  j2+1, k2+1) += Fx2 * Wy2       * Wz2;

      // jy
      if (D >= 1) yee.jy(i1,  j1,   k1)   += Fy1 * (1.0-Wx1) * (1.0-Wz1);
      if (D >= 2) yee.jy(i1+1,j1,   k1)   += Fy1 * Wx1       * (1.0-Wz1);
      if (D >= 3) yee.jy(i1  ,j1,   k1+1) += Fy1 * (1.0-Wx1) * Wz1;
      if (D >= 3) yee.jy(i1+1,j1,   k1+1) += Fy1 * Wx1       * Wz1;
      
      if (D >= 1) yee.jy(i2,  j2,   k2)   += Fy2 * (1.0-Wx2) * (1.0-Wz2);
      if (D >= 2) yee.jy(i2+1,j2,   k2)   += Fy2 * Wx2       * (1.0-Wz2);
      if (D >= 3) yee.jy(i2,  j2,   k2+1) += Fy2 * (1.0-Wx2) * Wz2;
      if (D >= 3) yee.jy(i2+1,j2,   k2+1) += Fy2 * Wx2       * Wz2;
                            
      // jz
      yee.jz(i1,  j1,   k1)   += Fz1 * (1.0-Wx1) * (1.0-Wy1);
      yee.jz(i1+1,j1,   k1)   += Fz1 * Wx1       * (1.0-Wy1);
      yee.jz(i1,  j1+1, k1)   += Fz1 * (1.0-Wx1) * Wy1;
      yee.jz(i1+1,j1+1, k1)   += Fz1 * Wx1       * Wy1;

      yee.jz(i2,  j2,   k2)   += Fz2 * (1.0-Wx2) * (1.0-Wy2);
      yee.jz(i2+1,j2,   k2)   += Fz2 * Wx2       * (1.0-Wy2);
      yee.jz(i2,  j2+1, k2)   += Fz2 * (1.0-Wx2) * Wy2;
      yee.jz(i2+1,j2+1, k2)   += Fz2 * Wx2       * Wy2;

    }

  }//end of loop over species

}



//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag<1,3>; // 1D3V
template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag<3,3>; // 3D3V

