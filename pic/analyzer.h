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


    // analysis lattice reference
    int ispc = 0;
    auto& analysis = cell.analysis[ispc];
    analysis.rho.clear();
    analysis.mgamma.clear();
    analysis.Vx.clear();
    analysis.Tx.clear();
    analysis.ekin.clear();

    // initialize pointers to particle arrays
    int nparts = cell.container.size();
      
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( cell.container.vel(i,0) );


    double gam;
    //double c = cell.cfl;
    double q = cell.container.qe; // TODO: split into species
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

      gam = sqrt(1.0 + 
          vel[0][n]*vel[0][n] + 
          vel[1][n]*vel[1][n] + 
          vel[2][n]*vel[2][n]);

      // --------------------------------------------------
      // general quantities
        
      yee.rho(i,j,k) += abs(q); // number density

      // --------------------------------------------------
      // particle-species quantities

      analysis.rho(i,j,k) += abs(q); // number density

      analysis.mgamma(i,j,k) += gam; // mean gamma

      analysis.Vx(i,j,k) += u0/gam; // bulk velocity

      // TODO Tx
      
      // kinetic energy
      // chi(u) = 1/2 m v.v
      analysis.ekin(i,j,k) += 0.5*abs(q)*u0*u0/gam/gam;


    }


    // normalize weight with number density
    for (size_t k=0; k<cell.NzMesh; k++)
    for (size_t j=0; j<cell.NyMesh; j++)
    for (size_t i=0; i<cell.NxMesh; i++) {
      analysis.mgamma(i,j,k) /= analysis.rho(i,j,k);
    }



    return;
  }


};



} // end of namespace pic
