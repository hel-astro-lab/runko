#pragma once

#include <algorithm>
#include <assert.h>

#include "cell.h"

//#include <fmt/format.h>
//#include <fmt/format.cc>
//#include <fmt/string.h>
//#include <fmt/ostream.h>


using std::min;
using std::max;


namespace pic {

class Depositer {
  public:


  void deposit( pic::PicCell& cell)
  {

    auto& yee = cell.getYee();
    yee.jx.clear();
    yee.jy.clear();
    yee.jz.clear();

    auto mins = cell.mins;

    // initialize pointers to particle arrays
    int nparts = cell.container.size();
      
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( cell.container.vel(i,0) );


    double invgam;
    double c = cell.cfl;
    double q = cell.container.qe;

    double x0, y0, z0, x1sp, x2sp, y1sp, y2sp, z1sp, z2sp;



    //int i1,i2,j1,j2,k1,k2;
    //double i1,i2,j1,j2,k1,k2;
    //double i1p1,i2p1,j1p1,j2p1,k1p1,k2p1;
    int i1,i2,j1,j2,k1,k2;
    int i1p1,i2p1,j1p1,j2p1,k1p1,k2p1;

    double xr, yr, zr;

    double Fx1, Fy1, Fz1, Fx2, Fy2, Fz2;
    double Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;
    double onemWx1, onemWy1, onemWz1, onemWx2, onemWy2, onemWz2;

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

      //q = weight*qe;
      //q = 1.0;

      // normalized location w.r.t. tile
      x1sp = x0         - mins[0];
      x2sp = loc[0][n]  - mins[0];
      y1sp = y0         - mins[1];
      y2sp = loc[1][n]  - mins[1];
      z1sp = z0         - mins[2];
      z2sp = loc[2][n]  - mins[2];


		  i1  = (int)floor( x1sp );
		  i2  = (int)floor( x2sp );
		  j1  = (int)floor( y1sp );
		  j2  = (int)floor( y2sp );
		  k1  = (int)floor( z1sp );
		  k2  = (int)floor( z2sp );

      // relay point; +1 is equal to +\Delta x
	    xr = min( (double)(min(i1,i2)+1), max( (double)max(i1,i2), 0.5*(x1sp+x2sp) ) );
	    yr = min( (double)(min(j1,j2)+1), max( (double)max(j1,j2), 0.5*(y1sp+y2sp) ) );
	    zr = min( (double)(min(k1,k2)+1), max( (double)max(k1,k2), 0.5*(z1sp+z2sp) ) );

      // twoD mode
      k1 = 0;
      k2 = 0;

	    // -q to include -j in the Ampere's equation
	    Fx1 = -q*(xr - x1sp);
	    Fy1 = -q*(yr - y1sp);
	    Fz1 = -q*(zr - z1sp);
	    
	    Wx1 = 0.5*(x1sp + xr) - i1;
	    Wy1 = 0.5*(y1sp + yr) - j1;
      //Wz1 = 0.5*(z1sp + zr) - k1;

	    Wx2 = 0.5*(x2sp + xr) - i2;
    	Wy2 = 0.5*(y2sp + yr) - j2;
		  //Wz2 = 0.5*(z2sp + zr) - k2;

      // twoD mode
      Wz1 = 0.0; 
		  Wz2 = 0.0;


	    Fx2 = -q*(x2sp-xr);
	    Fy2 = -q*(y2sp-yr);
	    Fz2 = -q*(z2sp-zr);

      onemWx1 = 1.0 - Wx1;
      onemWx2 = 1.0 - Wx2;
      onemWy1 = 1.0 - Wy1;
      onemWy2 = 1.0 - Wy2;
      onemWz1 = 1.0 - Wz1;
      onemWz2 = 1.0 - Wz2;

      i1p1 = i1 + 1;
      i2p1 = i2 + 1;
      j1p1 = j1 + 1;
      j2p1 = j2 + 1;
      //k1p1 = k1 + 1;
      //k2p1 = k2 + 1;

      //twoD mode
      k1p1 = 0;
      k2p1 = 0;


      //fmt::print("n={} i1:{} i2:{} Fx1:{} Fx2:{}\n",n,i1,i2,Fx1,Fx2);

      /*
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

      std::cout << " x1sp: " << x1sp;
      std::cout << " y1sp: " << y1sp;
      std::cout << " z1sp: " << z1sp;
      std::cout << " ||| ";
      std::cout << " x2sp: " << x2sp;
      std::cout << " y2sp: " << y2sp;
      std::cout << " z2sp: " << z2sp;
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
      */

      //--------------------------------------------------
      // index checking
      assert(i1   >= -3 && i1   < (int)cell.NxMesh+3);
      assert(j1   >= -3 && j1   < (int)cell.NyMesh+3);
      assert(k1   >= -3 && k1   < (int)cell.NzMesh+3);

      assert(i1p1 >= -3 && i1p1 < (int)cell.NxMesh+3);
      assert(j1p1 >= -3 && j1p1 < (int)cell.NyMesh+3);
      assert(k1p1 >= -3 && k1p1 < (int)cell.NzMesh+3);

      assert(i2   >= -3 && i2   < (int)cell.NxMesh+3);
      assert(j2   >= -3 && j2   < (int)cell.NyMesh+3);
      assert(k2   >= -3 && k2   < (int)cell.NzMesh+3);

      assert(i2p1 >= -3 && i2p1 < (int)cell.NxMesh+3);
      assert(j2p1 >= -3 && j2p1 < (int)cell.NyMesh+3);
      assert(k2p1 >= -3 && k2p1 < (int)cell.NzMesh+3);

      //--------------------------------------------------

      //yee.jx(i1,  j1,   k1)   += 1.0;
      //yee.jx(i1,  j1p1, k1)   += 1.0;
      //yee.jx(i2,  j2,   k2)   += 1.0;
      //yee.jx(i2,  j2p1, k2)   += 1.0;


      // jx
      yee.jx(i1,  j1,   k1)   += Fx1 * onemWy1 * onemWz1;
      yee.jx(i1,  j1p1, k1)   += Fx1 * Wy1     * onemWz1;
      yee.jx(i1,  j1,   k1p1) += Fx1 * onemWy1 * Wz1;
      yee.jx(i1,  j1p1, k1p1) += Fx1 * Wy1     * Wz1;

      yee.jx(i2,  j2,   k2)   += Fx2 * onemWy2 * onemWz2;
      yee.jx(i2,  j2p1, k2)   += Fx2 * Wy2     * onemWz2;
      yee.jx(i2,  j2,   k2p1) += Fx2 * onemWy2 * Wz2;
      yee.jx(i2,  j2p1, k2p1) += Fx2 * Wy2     * Wz2;

      // jy
      yee.jy(i1,  j1,   k1)   += Fy1 * onemWx1 * onemWz1;
      yee.jy(i1p1,j1,   k1)   += Fy1 * Wx1     * onemWz1;
      yee.jy(i1  ,j1,   k1p1) += Fy1 * onemWx1 * Wz1;
      yee.jy(i1p1,j1,   k1p1) += Fy1 * Wx1     * Wz1;
      
      yee.jy(i2,  j2,   k2)   += Fy2 * onemWx2 * onemWz2;
      yee.jy(i2p1,j2,   k2)   += Fy2 * Wx2     * onemWz2;
      yee.jy(i2,  j2,   k2p1) += Fy2 * onemWx2 * Wz2;
      yee.jy(i2p1,j2,   k2p1) += Fy2 * Wx2     * Wz2;
                            
      // jz
      yee.jz(i1,  j1,   k1)   += Fz1 * onemWx1 * onemWy1;
      yee.jz(i1p1,j1,   k1)   += Fz1 * Wx1     * onemWy1;
      yee.jz(i1,  j1p1, k1)   += Fz1 * onemWx1 * Wy1;
      yee.jz(i1p1,j1p1, k1)   += Fz1 * Wx1     * Wy1;

      yee.jz(i2,  j2,   k2)   += Fz2 * onemWx2 * onemWy2;
      yee.jz(i2p1,j2,   k2)   += Fz2 * Wx2     * onemWy2;
      yee.jz(i2,  j2p1, k2)   += Fz2 * onemWx2 * Wy2;
      yee.jz(i2p1,j2p1, k2)   += Fz2 * Wx2     * Wy2;

    }

  }



};


} // end of namespace pic
