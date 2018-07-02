#pragma once

#include <cmath> 
#include <assert.h>

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

    // cell limits
    auto mins = cell.mins;
    auto maxs = cell.maxs;

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

    // TODO: think SIMD (not possibly due to ijk writing to yee
    for(int n=n1; n<n2; n++) {

      x0 = loc[0][n];
      y0 = loc[1][n];
      z0 = loc[2][n];

      // grid coordinate location
      /*
      i = floor(x0);
      j = floor(y0);
      k = floor(z0);
      //k = 0; // TODO: explicit 2D dimensionality enforcement
      */

		  //i  = trunc( cell.NxMesh*(x0-mins[0])/(maxs[0]-mins[0]) );
		  //j  = trunc( cell.NyMesh*(y0-mins[1])/(maxs[1]-mins[1]) );
		  //k  = trunc( cell.NzMesh*(z0-mins[2])/(maxs[2]-mins[2]) );
        
      // fixed grid form assuming dx = 1
		  i  = floor( loc[0][n]-mins[0] );
		  j  = floor( loc[1][n]-mins[1] );
		  k  = floor( loc[2][n]-mins[2] );

      /*
      std::cout << "----------------------\n";
      std::cout << "cell ijk =( " << cell.my_i << "," << cell.my_j << ")\n";
      std::cout << "nx ny nz "    << cell.Nx << " " << cell.Ny << "\n";
      std::cout << "nxG nyG nzG " << cell.NxMesh << " " << cell.NyMesh << " " << cell.NzMesh << "\n";
      std::cout << "ijk =(" << i << "," << j << "," << k << ")\n";
      std::cout << "mins " << mins[0] << " " << mins[1] << " " << mins[2] << "\n";
      std::cout << "maxs " << maxs[0] << " " << maxs[1] << " " << maxs[2] << "\n";
      std::cout << "x "    << x0 << " " << y0 << " " << z0 << "\n";
      std::cout << " n = " << n << " " << n1 << " " << n2 << "\n";
      */

      assert(i >= 0 && i < cell.NxMesh);
      assert(j >= 0 && j < cell.NyMesh);
      assert(k >= 0 && k < cell.NzMesh);

      assert( x0 >= mins[0] && i < maxs[0] );
      assert( y0 >= mins[1] && j < maxs[1] );
      assert( z0 >= mins[2] && k < maxs[2] );

      u0 = vel[0][n];
      v0 = vel[1][n];
      w0 = vel[2][n];

      gam = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);

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
    for (size_t i=0; i<cell.NxMesh; i++)
    for (size_t j=0; j<cell.NyMesh; j++)
    for (size_t k=0; k<cell.NzMesh; k++)
      analysis.mgamma(i,j,k) /= analysis.rho(i,j,k);



    return;
  }


};



} // end of namespace pic
