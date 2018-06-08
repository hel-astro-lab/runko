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

    auto mins = cell.mins;
    auto maxs = cell.maxs;

    #pragma omp simd
    for(int n=n1; n<n2; n++) {

      // particle location in the grid
		  i  = trunc( cell.NxMesh*(loc[0][n]-mins[0])/(maxs[0]-mins[0]) );
		  dx = (loc[0][n]-mins[0]) - i;

		  j  = trunc( cell.NyMesh*(loc[1][n]-mins[1])/(maxs[1]-mins[1]) );
		  dy = (loc[1][n]-mins[1]) - j;

		  k  = trunc( cell.NzMesh*(loc[2][n]-mins[2])/(maxs[2]-mins[2]) );
		  dz = (loc[2][n]-mins[2]) - k;

      //std::cout << '\n';
      //std::cout << "x: " << loc[0][n] << " y: " << loc[1][n] << " z:" << loc[2][n] << '\n';
      //std::cout << " ijk " << i << "," << j << "," << k << '\n';
      //std::cout << " ds " << dx << "," << dy << "," << dz << '\n';
		  //l = i; // + iy*(j-1) + iz*(k-1);


      // TODO: 2D hack
      k = 0;
      dz = 0.0;
      int iz = 0;

      f = yee.ex(i,j,k) + yee.ex(i-1,j,k) +    dx*(yee.ex(i+1,j  ,k   ) - yee.ex(i-1,j  ,k  ));
      f+=                                      dy*(yee.ex(i  ,j+1,k   ) + yee.ex(i-1,j+1,k  )+
                                               dx*(yee.ex(i+1,j+1,k   ) - yee.ex(i-1,j+1,k  ))-f);
      g = yee.ex(i,j,k+1)+yee.ex(i-1,j,k+iz)+  dx*(yee.ex(i+1,j  ,k+iz) - yee.ex(i-1,j  ,k+iz));
      g+=                                      dy*(yee.ex(i  ,j+1,k+iz) + yee.ex(i-1,j+1,k+iz)
                                           +   dx*(yee.ex(i+1,j+1,k+iz) - yee.ex(i-1,j+1,k+iz))-g);
      ex[n] = (f+dz*(g-f))*0.5;

      f = yee.ey(i,j,k)+yee.ey(i,j-1,k)+       dy*(yee.ey(i  ,j+1,k   ) - yee.ey(i  ,j-1,k  ));
      f+=                                      dz*(yee.ey(i  ,j  ,k+iz) + yee.ey(i  ,j-1,k+iz)+      
                                               dy*(yee.ey(i  ,j+1,k+iz) - yee.ey(i  ,j-1,k+iz))-f);
      g = yee.ey(i+1,j,k)+yee.ey(i+1,j-1,k)+   dy*(yee.ey(i+1,j+1,k   ) - yee.ey(i+1,j-1,k   ));
      g+=                                      dz*(yee.ey(i+1,j  ,k+iz) + yee.ey(i+1,j-1,k+iz)+ 
                                               dy*(yee.ey(i+1,j+1,k+iz) - yee.ey(i+1,j-1,k+iz))-g);
      ey[n]=(f+dx*(g-f))*0.5;

      f = yee.ez(i,j,k)+yee.ez(i,j,k-iz)+      dz*(yee.ez(i  ,j  ,k+iz) - yee.ez(i  ,j  ,k-iz));
      f+=                                      dx*(yee.ez(i+1,j  ,k   ) + yee.ez(i+1,j  ,k-iz)+
                                               dz*(yee.ez(i+1,j  ,k+iz) - yee.ez(i+1,j  ,k-iz))-f);
      g = yee.ez(i,j+1,k)+ yee.ez(i,j+1,k-iz)+ dz*(yee.ez(i  ,j+1,k+iz) - yee.ez(i  ,j+1,k-iz));
      g+=                                      dx*(yee.ez(i+1,j+1,k   ) + yee.ez(i+1,j+1,k-iz)+
                                               dz*(yee.ez(i+1,j+1,k+iz) - yee.ez(i+1,j+1,k-iz))-g);
      ez[n]=(f+dy*(g-f))*0.5;

      f = yee.bx(i,j-1,k)  +yee.bx(i,j-1,k-iz )   +dz*(yee.bx(i,j-1,k+iz)   - yee.bx(i,j-1,k-iz));
      f = yee.bx(i,j,k)    +yee.bx(i,j,k-iz)      +dz*(yee.bx(i,j,k+iz)     - yee.bx(i,j,k-iz))+f+dy 
        *(yee.bx(i,j+1,k)  +yee.bx(i,j+1,k-iz)    +dz*(yee.bx(i,j+1,k+iz)   - yee.bx(i,j+1,k-iz))-f);
      g = yee.bx(i+1,j-1,k)+yee.bx(i+1,j-1,k-iz)  +dz*(yee.bx(i+1,j-1,k+iz) - yee.bx(i+1,j-1,k-iz));
      g = yee.bx(i+1,j,k)  +yee.bx(i+1,j,k-iz)    +dz*(yee.bx(i+1,j,k+iz)   - yee.bx(i+1,j,k-iz))
                                           +g     +dy*(yee.bx(i+1,j+1,k)    + yee.bx(i+1,j+1,k-iz)
                                                  +dz*(yee.bx(i+1,j+1,k+iz) - yee.bx(i+1,j+1,k-iz))-g);
      bx[n]=(f+dx*(g-f))*(.25*cinv);

      f = yee.by(i,j,k-iz)+yee.by(i-1,j,k-iz)     +dx*(yee.by(i+1,j,k-iz)   - yee.by(i-1,j,k-iz));
      f = yee.by(i,j,k)+yee.by(i-1,j,k)           +dx*(yee.by(i+1,j,k)      - yee.by(i-1,j,k))+f+dz 
        *(yee.by(i,j,k+iz)+yee.by(i-1,j,k+iz)     +dx*(yee.by(i+1,j,k+iz)   - yee.by(i-1,j,k+iz))-f);
      g = yee.by(i,j+1,k-iz)+yee.by(i-1,j+1,k-iz) +dx*(yee.by(i+1,j+1,k-iz) - yee.by(i-1,j+1,k-iz));
      g = yee.by(i,j+1,k)+yee.by(i-1,j+1,k)       +dx*(yee.by(i+1,j+1,k)    - yee.by(i-1,j+1,k))
                                               +g +dz*(yee.by(i,j+1,k+iz)   + yee.by(i-1,j+1,k+iz)
                                                  +dx*(yee.by(i+1,j+1,k+iz) - yee.by(i-1,j+1,k+iz))-g);
      by[n]=(f+dy*(g-f))*(.25*cinv);

      f = yee.bz(i-1,j,k)+yee.bz(i-1,j-1,k )      +dy*(yee.bz(i-1,j+1,k)    - yee.bz(i-1,j-1,k));
      f = yee.bz(i,j,k)+yee.bz(i,j-1,k)           +dy*(yee.bz(i,j+1,k)      - yee.bz(i,j-1,k))+f+dx 
        * (yee.bz(i+1,j,k)+yee.bz(i+1,j-1,k)      +dy*(yee.bz(i+1,j+1,k)    - yee.bz(i+1,j-1,k))-f);
      g = yee.bz(i-1,j, k+iz)+yee.bz(i-1,j-1,k+iz)+dy*(yee.bz(i-1,j+1,k+iz) - yee.bz(i-1,j-1,k+iz));
      g = yee.bz(i,j,k+iz)+yee.bz(i,j-1,k+iz )+    dy*(yee.bz(i,j+1,k+iz)   - yee.bz(i,j-1,k+iz))
                                                +g+dx*(yee.bz(i+1,j,k+iz)   + yee.bz(i+1,j-1,k+iz)
                                                  +dy*(yee.bz(i+1,j+1,k+iz) - yee.bz(i+1,j-1,k+iz))-g);
      bz[n]=(f+dz*(g-f))*(.25*cinv);
    }

    return;
  }


};




} // end of namespace pic
