#include "boris.h"

#include <cmath> 
#include "../../tools/signum.h"

//#include <fmt/format.h>
//#include <fmt/format.cc>
//#include <fmt/string.h>
//#include <fmt/ostream.h>

using toolbox::sign;


template<size_t D, size_t V>
void pic::BorisPusher<D,V>::solve(
    pic::Tile<D>& tile)
{

  // get reference to the Yee grid 
  //auto& yee = tile.get_yee();

  for(auto&& container : tile.containers) {
    int nparts = container.size();


    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );

    /*
    double* ex = &( (*tile.container.Epart)[0*nparts] );
    double* ey = &( (*tile.container.Epart)[1*nparts] );
    double* ez = &( (*tile.container.Epart)[2*nparts] );

    double* bx = &( (*tile.container.Bpart)[0*nparts] );
    double* by = &( (*tile.container.Bpart)[1*nparts] );
    double* bz = &( (*tile.container.Bpart)[2*nparts] );
    */

    double ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
    double bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

    double *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( container.Epart[0*nparts] );
    ey = &( container.Epart[1*nparts] );
    ez = &( container.Epart[2*nparts] );

    bx = &( container.Bpart[0*nparts] );
    by = &( container.Bpart[1*nparts] );
    bz = &( container.Bpart[2*nparts] );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    double u0, v0, w0;
    double u1, v1, w1;
    double g, f;


    double c = tile.cfl;
    double cinv = 1.0/c;

    // charge (sign only)
    double qm = sign(container.q);
    //double qm = container.q;
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
      for(size_t i=0; i<D; i++)
        loc[i][n] += vel[i][n]*g*c;


      //fmt::print("dx: {}\n",  vel[0][n]*g*c);
    }

  }//end of loop over species

}




//--------------------------------------------------
// explicit template instantiation

template class pic::BorisPusher<1,3>; // 1D3V
template class pic::BorisPusher<2,3>; // 2D3V
template class pic::BorisPusher<3,3>; // 3D3V

