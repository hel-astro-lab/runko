#pragma once

#include <cmath> 

#include "cell.h"
#include "../em-fields/fields.h"

#include "../tools/signum.h"

#include <iostream>


using toolbox::sign;


namespace pic {


/// General analyzator that computes moments for the particles inside the cells
//template<typename T>
class Analyzator {


  public:

  virtual void analyze( pic::PicCell& cell )
  {

    // Yee lattice reference
    auto& yee = cell.getYee();
    yee.rho.clear();
    yee.ekin.clear();
    //yee.jx1.clear();


    // initialize pointers to particle arrays
    int nparts = cell.container.size();
      
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( cell.container.vel(i,0) );


    double gam;
    double c = cell.cfl;
    double q = cell.container.qe;
    double x0, y0, z0;
    double u0, v0, w0;
    double i,j,k;


    // loop and check particles
    int n1 = 0;
    int n2 = nparts;

    #pragma omp simd 
    for(int n=n1; n<n2; n++) {

      x0 = loc[0][n];
      y0 = loc[1][n];
      z0 = loc[2][n];

      // grid coordinate location
      i = floor(x0);
      j = floor(y0);
      k = floor(z0);

      u0 = vel[0][n];
      v0 = vel[1][n];
      w0 = vel[2][n];

      gam = 1.0/sqrt(1.0 + 
          vel[0][n]*vel[0][n] + 
          vel[1][n]*vel[1][n] + 
          vel[2][n]*vel[2][n]);


      // number density
      yee.rho(i,j,k) += abs(q);

      // kinetic energy
      yee.ekin(i,j,k) += abs(q)*u0*u0;
      //std::cout << abs(q)*u0*u0 << "ijk:" << i << " " << j << " " << k << '\n';
    }


    return;
  }


};



} // end of namespace pic
