#pragma once

#include <cmath> 

//#include <fmt/format.h>
//#include <fmt/format.cc>
//#include <fmt/string.h>
//#include <fmt/ostream.h>

#include "../tools/signum.h"
#include "../definitions.h"

#include "tile.h"

namespace pic {

  using toolbox::sign;

class Pusher
{

  public:

  /*! \brief Push particles
   *
   */
  void solve(pic::Tile<2>& tile)
  {

    // get reference to the Yee grid 
    //auto& yee = tile.getYee();

    for (size_t ispc=0; ispc<tile.Nspecies(); ispc++) {
      ParticleBlock& container = tile.get_container(ispc);

      int nparts = container.size();


      // initialize pointers to particle arrays
      Realf* loc[3];
      for( int i=0; i<3; i++)
        loc[i] = &( container.loc(i,0) );

      Realf* vel[3];
      for( int i=0; i<3; i++)
        vel[i] = &( container.vel(i,0) );

      /*
      Realf* ex = &( (*tile.container.Epart)[0*nparts] );
      Realf* ey = &( (*tile.container.Epart)[1*nparts] );
      Realf* ez = &( (*tile.container.Epart)[2*nparts] );

      Realf* bx = &( (*tile.container.Bpart)[0*nparts] );
      Realf* by = &( (*tile.container.Bpart)[1*nparts] );
      Realf* bz = &( (*tile.container.Bpart)[2*nparts] );
      */

      Realf ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
      Realf bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

      Realf *ex, *ey, *ez, *bx, *by, *bz;
      ex = &( container.Epart[0][0] );
      ey = &( container.Epart[1][0] );
      ez = &( container.Epart[2][0] );

      bx = &( container.Bpart[0][0] );
      by = &( container.Bpart[1][0] );
      bz = &( container.Bpart[2][0] );

      // loop over particles
      int n1 = 0;
      int n2 = nparts;

      Realf u0, v0, w0;
      Realf u1, v1, w1;
      Realf g, f;


      Realf c = tile.cfl;
      Realf cinv = 1.0/c;

      // charge (sign only)
      Realf qm = sign(container.q);
      //Realf qm = container.q;
      //std::cout << " qm = " << qm << " ispc: " << ispc << '\n';


      //TODO: SIMD
      //TODO: vectorize with matrix multiplication form & Eigen
      for(int n=n1; n<n2; n++) {

        //--------------------------------------------------
        // Boris algorithm

        // read particle-specific fields
        ex0 = ex[n]*(0.5*qm);
        ey0 = ey[n]*(0.5*qm);
        ez0 = ez[n]*(0.5*qm);

        bx0 = bx[n]*(0.5*qm*cinv);
        by0 = by[n]*(0.5*qm*cinv);
        bz0 = bz[n]*(0.5*qm*cinv);


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
	  		u1 = (u0 + v0*bz0 - w0*by0)*f;
	  		v1 = (v0 + w0*bx0 - u0*bz0)*f;
	  		w1 = (w0 + u0*by0 - v0*bx0)*f;

        // second half of magnetic rotation & electric acceleration
	  		u0 = u0 + v1*bz0 - w1*by0 + ex0;
	  		v0 = v0 + w1*bx0 - u1*bz0 + ey0;
	  		w0 = w0 + u1*by0 - v1*bx0 + ez0;

        // normalized 4-velocity advance
        vel[0][n] = u0*cinv;
        vel[1][n] = v0*cinv;
        vel[2][n] = w0*cinv;

        /*
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        std::cout << "n=" << n;
        std::cout << "\n";

        std::cout << " vx: " <<  vel[0][n];
        std::cout << " vy: " <<  vel[1][n];
        std::cout << " vz: " <<  vel[2][n];
        std::cout << "\n";

        std::cout << " ex: " <<  ex[n];
        std::cout << " ey: " <<  ey[n];
        std::cout << " ez: " <<  ez[n];
        std::cout << " bx: " <<  bx[n];
        std::cout << " by: " <<  by[n];
        std::cout << " bz: " <<  bz[n];
        std::cout << "\n";
        */


        // position advance
	  	  g = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
        //TODO: note the explicit 2D dimensionality enforcement
        //for(int i=0; i<2; i++)
        //  loc[i][n] += vel[i][n]*g*c;

        loc[0][n] += vel[0][n]*g*c;
        loc[1][n] += vel[1][n]*g*c;
        //loc[2][n] += vel[2][n]*g*c;

        //fmt::print("dx: {}\n",  vel[0][n]*g*c);
      }

    }//end of loop over species
  }

};





} // end of namespace pic
