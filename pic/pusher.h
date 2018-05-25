#pragma once

#include <cmath> 

//#include <fmt/format.h>
//#include <fmt/format.cc>
//#include <fmt/string.h>
//#include <fmt/ostream.h>


#include "cell.h"

namespace pic {


class Pusher
{

  public:

  /*! \brief Push particles
   *
   */
  void solve(pic::PicCell& cell)
  {

    // get reference to the Yee grid 
    //auto& yee = cell.getYee();

    int nparts = cell.container.size();


    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( cell.container.vel(i,0) );

    /*
    double* ex = &( (*cell.container.Epart)[0*nparts] );
    double* ey = &( (*cell.container.Epart)[1*nparts] );
    double* ez = &( (*cell.container.Epart)[2*nparts] );

    double* bx = &( (*cell.container.Bpart)[0*nparts] );
    double* by = &( (*cell.container.Bpart)[1*nparts] );
    double* bz = &( (*cell.container.Bpart)[2*nparts] );
    */

    double ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
    double bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

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

    double u0, v0, w0;
    double u1, v1, w1;
    double g, f;


    double c = cell.cfl;
    double cinv = 1.0/c;

    // electrons
    double qm = -1.0;



    #pragma omp simd
    for(int n=n1; n<n2; n++) {

      // read particle-specific fields
      ex0 = ex[n]*(0.5*qm);
      bx0 = bx[n]*(0.5*qm*cinv);

      // first half electric acceleration
      u0 = c*vel[0][n] + ex0;
      v0 = c*vel[1][n] + ey0;
      w0 = c*vel[2][n] + ez0;

      // first half magnetic rotation
      g = c/sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
      bx0 *= g;
      by0 *= g;
      bz0 *= g;

      f = 2.0/(1.0 + bx0*bx0 + by0*by0 + bz0*bz0);
			u1 = (u0 + v0*bz0 -w0*by0)*f;
			v1 = (v0 + w0*bx0 -u0*bz0)*f;
			w1 = (w0 + u0*by0 -v0*bx0)*f;

      // second half of magnetic rotation & electric acceleration
			u0 = u0 + v1*bz0 - w1*by0 + ex0;
			v0 = v0 + w1*bx0 - u1*bz0 + ey0;
			w0 = w0 + u1*by0 - v1*bx0 + ez0;

      // normalized 4-velocity advance
      vel[0][n] = u0*cinv;
      vel[1][n] = v0*cinv;
      vel[2][n] = w0*cinv;


      // position advance
		  //g = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
      g = 1.0;
      for(int i=0; i<3; i++)
        loc[i][n] += vel[i][n]*g*c;


      //fmt::print("dx: {}\n",  vel[0][n]*g*c);
    }

  }


};





} // end of namespace pic
