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

  virtual void analyze( pic::Tile<2>& tile )
  {

    // Yee lattice reference
    auto& yee = tile.getYee();
    yee.rho.clear();

    // tile limits
    auto mins = tile.mins;
    auto maxs = tile.maxs;

    // analysis lattice reference
    for (size_t ispc=0; ispc<tile.Nspecies(); ispc++) {
      ParticleBlock& container = tile.get_container(ispc);
      auto& analysis = tile.analysis[ispc];

      analysis.rho.clear();
      analysis.mgamma.clear();
      analysis.Vx.clear();
      analysis.Tx.clear();
      analysis.ekin.clear();


      // initialize pointers to particle arrays
      int nparts = container.size();
        
      Realf* loc[3];
      for( int i=0; i<3; i++)
        loc[i] = &( container.loc(i,0) );

      Realf* vel[3];
      for( int i=0; i<3; i++)
        vel[i] = &( container.vel(i,0) );


      Realf gam;
      //Realf c = tile.cfl;
      Realf q = container.q; // TODO: split into species
      Realf x0, y0, z0;
      Realf u0, v0, w0;
      //Realf i,j,k;
      int i,j,k;


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

		    //i  = trunc( tile.NxMesh*(x0-mins[0])/(maxs[0]-mins[0]) );
		    //j  = trunc( tile.NyMesh*(y0-mins[1])/(maxs[1]-mins[1]) );
		    //k  = trunc( tile.NzMesh*(z0-mins[2])/(maxs[2]-mins[2]) );
          
        // fixed grid form assuming dx = 1
		    i  = (int)floor( loc[0][n]-mins[0] );
		    j  = (int)floor( loc[1][n]-mins[1] );
		    //k  = (int)floor( loc[2][n]-mins[2] );
        k = 0; // 2D hack

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

        assert(i >= 0 && i < static_cast<int>(tile.mesh_lengths[0]) );
        assert(j >= 0 && j < static_cast<int>(tile.mesh_lengths[1]) );
        //assert(k >= 0 && k < static_cast<int>(tile.mesh_lengths[2]) );

        assert( x0 >= mins[0] && x0 < maxs[0] );
        assert( y0 >= mins[1] && y0 < maxs[1] );
        //assert( z0 >= mins[2] && z0 < maxs[2] );

        u0 = vel[0][n];
        v0 = vel[1][n];
        w0 = vel[2][n];

        gam = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);

        // --------------------------------------------------
        // general quantities
        yee.rho(i,j,k) += abs(q); // total number density

        analysis.rho(i,j,k) += abs(q); // number density per species
          

        // --------------------------------------------------
        // particle-species quantities
        analysis.mgamma(i,j,k) += gam; // mean gamma

        analysis.Vx(i,j,k) += u0/gam; // bulk velocity

        // TODO Tx
        
        // kinetic energy
        // chi(u) = 1/2 m v.v
        analysis.ekin(i,j,k) += 0.5*abs(q)*(u0*u0 + v0*v0 + w0*w0)/gam/gam;

      }

      // normalize weight with number density
      for (size_t i=0; i<tile.mesh_lengths[0]; i++)
      for (size_t j=0; j<tile.mesh_lengths[1]; j++)
        analysis.mgamma(i,j,0) /= analysis.rho(i,j,0);

    }


    return;
  }


};



} // end of namespace pic
