#pragma once

#include <cmath> 


#include "cell.h"

namespace pic {

class ParticleFieldInterpolator 
{

  public:

  /*! \brief interpolate electromagnetic fields to particle locations
   *
   */
  void solve(pic::PicCell& cell)
  {

    // get reference to the Yee grid 
    auto& yee = cell.getYee();

    int nparts = cell.container.size();


    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    // 1-d arrays
    //double* ex = &( (*cell.container.Epart)[0*nparts] );
    //double* ey = &( (*cell.container.Epart)[1*nparts] );
    //double* ez = &( (*cell.container.Epart)[2*nparts] );

    // multiD array version
    //double *efield[3], *bfield[3];
    //for( int i=0; i<3; i++) {
    //  efield[i] = &( cell.container.Epart[i][0] );
    //  bfield[i] = &( cell.container.Bpart[i][0] );
    //}
      
    double *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( cell.container.Epart[0][0] );
    ey = &( cell.container.Epart[1][0] );
    ez = &( cell.container.Epart[2][0] );
    bx = &( cell.container.Bpart[0][0] );
    by = &( cell.container.Bpart[1][0] );
    bz = &( cell.container.Bpart[2][0] );


    // loop over particles
    int n1 = 0;
    int n2 = nparts;


    double c = 1.0; // cfl
    double cinv = 1.0/c;
    //double qm = 2.0;


    int i=0, j=0, k=0;
    //int ix=0, iy=0, iz=0;
    double dx=0.0, dy=0.0, dz=0.0;
    double f,g;

    double xmin = 0.0; // block starting location
    double lenx = 1.0; // grid cell size in x dir

    #pragma omp simd
    for(int n=n1; n<n2; n++) {

      // particle location in the grid
		  i  = trunc( (loc[0][n]-xmin)/lenx );
		  dx = (loc[0][n]-xmin)*lenx - i;

		  j  = trunc( loc[1][n] );
		  dy = loc[1][n] - j;

		  k  = trunc( loc[2][n] );
		  dz = loc[2][n] - k;

		  //l = i; // + iy*(j-1) + iz*(k-1);

      // ex component
  		//f = ex[l] + ex[l-ix] + dx*(ex[l+ix] - ex[l-ix]);
		  //f = f + dy*(ex[l+iy] + ex[l-ix+iy] + dx*(ex[l+ix+iy] - ex[l-ix+iy]) - f);
  		//g = ex[l+iz] + ex[l-ix+iz] + dx*(ex[l+ix+iz] - ex[l-ix+iz]);
  		//g = g + dy*(ex[l+iy+iz] + ex[l-ix+iy+iz] + 
      //        dx*(ex[l+ix+iy+iz] - ex[l-ix+iy+iz]) - g);
		  //ex0 = (f + dz*(g - f))*(.25*qm);

  		f = yee.ex(i,j,k) + yee.ex(i-1,j,k) + dx*(yee.ex(i+1,j,k) - yee.ex(i-1,j,k));
      f += dy*( yee.ex(i,j+1,k) + yee.ex(i-1,j+1,k) + dx*(yee.ex(i+1, j+1, k) - yee.ex(i-1, j+1, k) ) - f);
      g = yee.ex(i,j,k+1) + yee.ex(i-1,j,k+1) + dx*(yee.ex(i+1,j,k+1) - yee.ex(i-1,j,k+1));
      g += dy*(yee.ex(i,j+1,k+1) + yee.ex(i-1,j+1,k+1) 
         + dx*(yee.ex(i+1,j+1,k+1) - yee.ex(i-1,j+1,k+1)) - g);
      ex[n] = (f + dz*(g-f))*0.5;

      //ex[n] = f*(0.25*qm); // normalize in pusher
      //ex[n] = f*0.5;


      // rest of the components TODO
      ey[n] = 0.0;
      ez[n] = 0.0;

      // bx component
		  // f=bx(l-iy,1,1)+bx(l-iy-iz,1,1)+dz*(bx(l-iy+iz,1,1)-bx(l-iy-iz,1,1))
		  // f=bx(l,1,1)+bx(l-iz,1,1)+dz*(bx(l+iz,1,1)-bx(l-iz,1,1))+f+dy &
		  // * (bx(l+iy,1,1)+bx(l+iy-iz,1,1)+dz*(bx(l+iy+iz,1,1)-bx(l+iy &
		  // -iz,1,1))-f)
		  // g=bx(l+ix-iy,1,1)+bx(l+ix-iy-iz,1,1)+dz*(bx(l+ix-iy+iz,1,1) &
		  // -bx(l+ix-iy-iz,1,1))
		  // g=bx(l+ix,1,1)+bx(l+ix-iz,1,1)+dz*(bx(l+ix+iz,1,1)-bx(l+ix-iz,1 &
		  // ,1))+g+dy*(bx(l+ix+iy,1,1)+bx(l+ix+iy-iz,1,1)+dz*(bx(l+ix &
		  // +iy+iz,1,1)-bx(l+ix+iy-iz,1,1))-g)
		  // 
		  // bx0=(f+dx*(g-f))*(.125*qm*cinv)

      f = yee.bx(i,j-1,k) + yee.bx(i,j-1,k-1) + dz*(yee.bx(i,j-1,k+1) - yee.bx(i,j-1,k-1));
      f = yee.bx(i,j,k) + yee.bx(i,j,k-1) + dz*(yee.bx(i,j,k+1) - yee.bx(i,j,k-1)) + f
          + dy*(yee.bx(i,j+1,k) + yee.bx(i,j+1,k-1) + dz*(yee.bx(i,j+1,k+1) - yee.bx(i,j+1,k-1)) -f);
      g = yee.bx(i+1,j-1,k) + yee.bx(i+1,j-1,k-1) + dz*(yee.bx(i+1,j-1,k+1) - yee.bx(i+1,j-1,k-1));
      g = yee.bx(i+1,j,k) + yee.bx(i+1,j,k-1) + dz*(yee.bx(i+1,j,k+1) - yee.bx(i+1,j,k-1)) + g + 
        dy*(yee.bx(i+1,j+1,k) + yee.bx(i+1,j+1,k-1) + dz*(yee.bx(i+1,j+1,k+1)-yee.bx(i+1,j+1,k-1))-g);
      bx[n] = (f + dx*(g-f))*(0.25*cinv);


      // rest of the components TODO
      by[n] = 0.0;
      bz[n] = 0.0;

    }

    return;
  }


};






} // end of namespace pic
