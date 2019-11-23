#pragma once

#include <cmath> 
#include <cassert>

#include "tile.h"
#include "../em-fields/tile.h"

#include "../tools/signum.h"

#include <iostream>


using toolbox::sign;


namespace pic {


/// General analyzator that computes moments for the particles inside the tiles
//template<typename T>
class Analyzator {


  public:

  Analyzator() {};

  virtual ~Analyzator() = default;

  template<int D>
  void analyze( pic::Tile<D>& tile )
  {

    // Yee lattice reference
    auto& yee = tile.get_yee();
    yee.rho.clear();

    // tile limits
    auto mins = tile.mins;
    auto maxs = tile.maxs;

    // analysis lattice reference
    for (int ispc=0; ispc<tile.Nspecies(); ispc++) {
      ParticleContainer& container = tile.get_container(ispc);
      auto& analysis = tile.analysis[ispc];

      analysis.clear();


      // initialize pointers to particle arrays
      int nparts = container.size();
        
      real_prtcl* loc[3];
      for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

      real_prtcl* vel[3];
      for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );


      //double c = tile.cfl;
      real_long gam;
      real_long q = container.q; // TODO: split into species
      real_long x0, y0, z0;
      real_long u0, v0, w0;
      int i,j,k;


      // loop and check particles
      int n1 = 0;
      int n2 = nparts;

      // TODO: think SIMD (not possible atm due to ijk writing to yee
      for(int n=n1; n<n2; n++) {
          
        // grid coordinate location; cast to double for the duration of this algorithm
        x0 = static_cast<real_long>( loc[0][n] );
        y0 = static_cast<real_long>( loc[1][n] );
        z0 = static_cast<real_long>( loc[2][n] );
          
        // fixed grid form assuming dx = 1
  	    i = D >= 1 ? static_cast<int>(floor( x0 - mins[0] ) ) : 0;
  	    j = D >= 2 ? static_cast<int>(floor( y0 - mins[1] ) ) : 0;
  	    k = D >= 3 ? static_cast<int>(floor( z0 - mins[2] ) ) : 0;

        /*
        std::cout << "----------------------\n";
        std::cout << "tile ijk =( " << tile.my_i << "," << tile.my_j << ")\n";
        std::cout << "nx ny nz "    << tile.Nx << " " << tile.Ny << "\n";
        std::cout << "nxG nyG nzG " << tile.NxMesh << " " << tile.NyMesh << " " << tile.NzMesh << "\n";
        std::cout << "ijk =(" << i << "," << j << "," << k << ")\n";
        std::cout << "mins " << mins[0] << " " << mins[1] << " " << mins[2] << "\n";
        std::cout << "maxs " << maxs[0] << " " << maxs[1] << " " << maxs[2] << "\n";
        std::cout << "x "    << x0 << " " << y0 << " " << z0 << "\n";
        std::cout << " n = " << n << " " << n1 << " " << n2 << "\n";
        */

        bool debug_flag = false;
        if(D >= 1){ if(! (i >= 0 && i < static_cast<int>(tile.mesh_lengths[0]))) debug_flag = true;}
        if(D >= 2){ if(! (j >= 0 && j < static_cast<int>(tile.mesh_lengths[1]))) debug_flag = true;}
        if(D >= 3){ if(! (k >= 0 && k < static_cast<int>(tile.mesh_lengths[2]))) debug_flag = true;}

        if(D >= 1){ if(! (x0 >= mins[0] && x0 < maxs[0]) ) debug_flag = true;}
        if(D >= 2){ if(! (y0 >= mins[1] && y0 < maxs[1]) ) debug_flag = true;}
        if(D >= 3){ if(! (z0 >= mins[2] && z0 < maxs[2]) ) debug_flag = true;}


        if (debug_flag) {

          std::cout << "--------------------------------------------------\n";
          std::cout << "n=" << n;
          std::cout << " i: " << i;
          std::cout << " j: " << j;
          std::cout << " k: " << k;
          std::cout << "\n";

          std::cout << " mins0: " << mins[0];
          std::cout << " mins1: " << mins[1];
          std::cout << " mins2: " << mins[2];

          std::cout << " maxs0: " << maxs[0];
          std::cout << " maxs1: " << maxs[1];
          std::cout << " maxs2: " << maxs[2];

          std::cout << " x0: " << x0;
          std::cout << " y0: " << y0;
          std::cout << " z0: " << z0;
          std::cout << " u0: " << vel[0][n];
          std::cout << " v0: " << vel[1][n];
          std::cout << " w0: " << vel[2][n];
          std::cout << "\n";
          
          std::cout << std::flush;
          // always fail
          assert(false);
        }

        u0 = static_cast<real_long>(vel[0][n]);
        v0 = static_cast<real_long>(vel[1][n]);
        w0 = static_cast<real_long>(vel[2][n]);

        gam = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);

        // --------------------------------------------------
        // general quantities
        real_long mass = abs(q);
        yee.rho(i,j,k) += mass; // total number density

        // --------------------------------------------------
        // particle-species quantities
        analysis.rho(i,j,k)     += mass; // number density per species

        analysis.edens(i,j,k)   += gam*mass; // energy density

        // momentum density
        analysis.momx(i,j,k)    += u0*mass;
        analysis.momy(i,j,k)    += v0*mass;
        analysis.momz(i,j,k)    += w0*mass;

        // pressure density
        analysis.pressx(i,j,k)  += u0*u0*mass/gam;
        analysis.pressy(i,j,k)  += v0*v0*mass/gam;
        analysis.pressz(i,j,k)  += w0*w0*mass/gam;

        // off-diagonal shear terms
        analysis.shearxy(i,j,k) += u0*v0*mass/gam;
        analysis.shearxz(i,j,k) += u0*w0*mass/gam;
        analysis.shearyz(i,j,k) += v0*w0*mass/gam;

        analysis.Vx(i,j,k)      += u0/gam; // bulk velocity x
        analysis.Vy(i,j,k)      += v0/gam; // bulk velocity y
        analysis.Vz(i,j,k)      += w0/gam; // bulk velocity z

        analysis.temp(i,j,k)    += gam - 1.0; // temperature XXX is it?

      }

      // normalize weight with number density
      //for (size_t i=0; i<tile.mesh_lengths[0]; i++)
      //for (size_t j=0; j<tile.mesh_lengths[1]; j++)
      //  analysis.mgamma(i,j,0) /= analysis.rho(i,j,0);

    }


    return;
  }


};



} // end of namespace pic
