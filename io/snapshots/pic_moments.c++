#include <cmath>

#include "io/snapshots/pic_moments.h"
#include "external/ezh5/src/ezh5.hpp"
#include "pic/particle.h"
#include "pic/tile.h"
#include "tools/signum.h"
#include "tools/limit.h"

// TODO turning compiler warnings off temporarily in this file since 
//      error printing in debug mode accesses mins/maxs outside boundaries
// NOTE remember to remove pragma pop at the end of the file when done with these.
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Warray-bounds"

using ezh5::File;


template<size_t D>
inline void h5io::PicMomentsWriter<D>::read_tiles(
    corgi::Grid<D>& grid)
{
  // clear target arrays
  for(auto& arr : arrs ) arr.clear();

  // target arrays
  auto& dense  = arrs[0];
  auto& densp =  arrs[1];

  auto& Vxe    = arrs[2];
  auto& Vye    = arrs[3];
  auto& Vze    = arrs[4];

  auto& Vxp    = arrs[5];
  auto& Vyp    = arrs[6];
  auto& Vzp    = arrs[7];

  auto& pressx = arrs[8];
  auto& pressy = arrs[9];
  auto& pressz = arrs[10];

  auto& shearxy= arrs[11];
  auto& shearxz= arrs[12];
  auto& shearyz= arrs[13];

  auto& densx =  arrs[14];

  // local variables
  double gam, xene, mass, wgt;
  double x0, y0, z0;
  double u0, v0, w0;
  int nparts;
  int i,j,k;
  int iff,jff,kff;


  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    auto mins = tile.mins;
    auto maxs = tile.maxs;

    // update also gs
    auto& gs = tile.get_grids();
    gs.rho.clear();

    // loop over species
    for (int ispc=0; ispc<tile.Nspecies(); ispc++) {
      auto& container = tile.get_container(ispc);
      mass = container.m; // species mass
      nparts = container.size();

      if(nparts <= 0) continue; // skip zero containers

      float* loc[3];
      for( i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

      float* vel[3];
      for( i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

      float* ch;
      ch = &( container.wgt(0) );


      // loop and search over all particles
      int n1 = 0;
      int n2 = nparts;
      for(int n=n1; n<n2; n++) {

        // prtcl coordinate location; cast to double for the duration of this algorithm
        x0  = static_cast<double>( loc[0][n] );
        y0  = static_cast<double>( loc[1][n] );
        z0  = static_cast<double>( loc[2][n] );
        wgt = static_cast<double>(     ch[n] );

        u0  = static_cast<double>( vel[0][n] );
        v0  = static_cast<double>( vel[1][n] );
        w0  = static_cast<double>( vel[2][n] );

        gam = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);
        xene = sqrt(      u0*u0 + v0*v0 + w0*w0);

#ifdef DEBUG
        // capture NaNs
        assert(!std::isnan(x0));
        assert(!std::isnan(y0));
        assert(!std::isnan(z0));
        assert(!std::isnan(wgt));

        // check that particles are communicated properly
        bool flagx = D>=1 ? (x0-mins[0] >= -3.0) && (x0 <= maxs[0] +2.0 ) : true;
        bool flagy = D>=2 ? (y0-mins[1] >= -3.0) && (y0 <= maxs[1] +2.0 ) : true;
        bool flagz = D>=3 ? (z0-mins[2] >= -3.0) && (z0 <= maxs[2] +2.0 ) : true;

        if( !flagx || !flagy || !flagz) {
          std::cout << "ERR IN MOM:" << std::endl;
          std::cerr << " minx: " << mins[0];
          std::cerr << " miny: " << mins[1];
          std::cerr << " minz: " << mins[2] << std::endl;
          std::cerr << " maxx: " << maxs[0];
          std::cerr << " maxy: " << maxs[1];
          std::cerr << " maxz: " << maxs[2] << std::endl;
          std::cerr << " x: " << x0-mins[0];
          std::cerr << " y: " << y0-mins[1];
          std::cerr << " z: " << z0-mins[2] << std::endl;
          std::cerr << " vx: " << u0;
          std::cerr << " vy: " << v0;
          std::cerr << " vz: " << w0 << std::endl;
          std::cerr << " fx: " << flagx;
          std::cerr << " fy: " << flagy;
          std::cerr << " fz: " << flagz;
          assert(false);
        }
#endif

        // rel prtcl index; assuming dx = 1; tile coordinates
        // limit to 0 Nx-1 just in case to avoid crashes
        iff = D >= 1 ? limit( floor(x0-mins[0]), -3., maxs[0]-mins[0] +2.) : 0;
        jff = D >= 2 ? limit( floor(y0-mins[1]), -3., maxs[1]-mins[1] +2.) : 0;
        kff = D >= 3 ? limit( floor(z0-mins[2]), -3., maxs[2]-mins[2] +2.) : 0;

        // update rho arrays; this is interpreted as mass density
        gs.rho(iff,jff,kff) += mass*wgt;

        //-------------------------------------------------- 

        // full prtcl index; assuming dx = 1; global grid coordinates
        // reduce by a factor of stride 
        i = D >= 1 ? limit( floor(x0/stride), 0.0, double(nx)-1.0) : 0;
        j = D >= 2 ? limit( floor(y0/stride), 0.0, double(ny)-1.0) : 0;
        k = D >= 3 ? limit( floor(z0/stride), 0.0, double(nz)-1.0) : 0;


        //--------------------------------------------------
        // next, physics quantities

        // number density 
        if(ispc == 0) dense(i,j,k) += wgt; 
        if(ispc == 1) densp(i,j,k) += wgt; 
        //if(ispc == 2) densx(i,j,k) += wgt; 
        if(ispc == 2) densx(i,j,k) += wgt*xene;  // energy density for photons

        // bulk flows
        if(ispc == 0) {
          Vxe(i,j,k) += wgt*u0/gam; 
          Vye(i,j,k) += wgt*v0/gam;
          Vze(i,j,k) += wgt*w0/gam;
        }
        if(ispc == 1) {
          Vxp(i,j,k) += wgt*u0/gam; 
          Vyp(i,j,k) += wgt*v0/gam;
          Vzp(i,j,k) += wgt*w0/gam;
        }

        // momentum density (flux of mass)
        // momx(i,j,k)    += u0*mass;
        // momy(i,j,k)    += v0*mass;
        // momz(i,j,k)    += w0*mass;
          
        // pressure (flux of momentum)
        pressx(i,j,k)  += wgt*u0*u0*mass/gam;
        pressy(i,j,k)  += wgt*v0*v0*mass/gam;
        pressz(i,j,k)  += wgt*w0*w0*mass/gam;

        // off-diagonal shear terms 
        shearxy(i,j,k) += wgt*u0*v0*mass/gam;
        shearxz(i,j,k) += wgt*u0*w0*mass/gam;
        shearyz(i,j,k) += wgt*v0*w0*mass/gam;

      } // end of prtcls
    } // end of species
  } // end of tiles

  // normalize bulk flow with local number of particles
  for(int k=0; k<nz; k++)
  for(int j=0; j<ny; j++)
  for(int i=0; i<nx; i++) {
    if(dense(i,j,k) > 0) {
      Vxe(i,j,k) /= dense(i,j,k);
      Vye(i,j,k) /= dense(i,j,k);
      Vze(i,j,k) /= dense(i,j,k);
    }

    if(densp(i,j,k) > 0) {
      Vxp(i,j,k) /= densp(i,j,k);
      Vyp(i,j,k) /= densp(i,j,k);
      Vzp(i,j,k) /= densp(i,j,k);
    }
  }
}



template<size_t D>
inline bool h5io::PicMomentsWriter<D>::write(
    corgi::Grid<D>& grid, int lap)
{
  read_tiles(grid);
  mpi_reduce_snapshots(grid);

  if( grid.comm.rank() == 0 ) {

    // build filename
    std::string full_filename = 
      fname + "/" +
      file_name + 
      "_" +
      std::to_string(lap) +
      extension;
    //std::cout << "QW: " << full_filename << std::endl;

    // open file and write
    File file(full_filename, H5F_ACC_TRUNC);
    file["Nx"] = arrs[0].Nx;
    file["Ny"] = arrs[0].Ny;
    file["Nz"] = arrs[0].Nz;


    file["dense"] = arrs[0].serialize();
    file["densp"] = arrs[1].serialize();
    file["densx"] = arrs[14].serialize(); // NOTE index

    file["Vxe"]   = arrs[2].serialize();
    file["Vye"]   = arrs[3].serialize();
    file["Vze"]   = arrs[4].serialize();

    file["Vxp"]   = arrs[5].serialize();
    file["Vyp"]   = arrs[6].serialize();
    file["Vzp"]   = arrs[7].serialize();

    file["pressx"]   = arrs[8].serialize();
    file["pressy"]   = arrs[9].serialize();
    file["pressz"]   = arrs[10].serialize();

    file["shearxy"]   = arrs[11].serialize();
    file["shearxz"]   = arrs[12].serialize();
    file["shearyz"]   = arrs[13].serialize();

  }

  return true;
}




//--------------------------------------------------
// explicit template class instantiations
template class h5io::PicMomentsWriter<1>;
template class h5io::PicMomentsWriter<2>;
template class h5io::PicMomentsWriter<3>;

#pragma GCC diagnostic pop
