#include <cmath>

#include "pic_moments.h"
#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../pic/particle.h"
#include "../../pic/tile.h"
#include "../../tools/signum.h"
#include "../../tools/limit.h"


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


  // local variables
  real_long gam;
  real_long mass;
  real_long charge;
  real_long x0, y0, z0;
  real_long u0, v0, w0;
  int nparts;
  int i,j,k;
  int iff,jff,kff;


  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
    auto mins = tile.mins;

    // update also yee
    auto& yee = tile.get_yee();
    yee.rho.clear();


    // loop over species
    for (int ispc=0; ispc<tile.Nspecies(); ispc++) {
      auto& container = tile.get_container(ispc);
      charge = container.q; // species charge
      nparts = container.size();

      real_prtcl* loc[3];
      for( i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

      real_prtcl* vel[3];
      for( i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

      //real_prtcl* ch;
      //ch = &( container.wgt(0) );


      // loop and search over all particles
      int n1 = 0;
      int n2 = nparts;
      for(int n=n1; n<n2; n++) {
        mass = std::abs(charge);

        // prtcl coordinate location; cast to double for the duration of this algorithm
        x0 = static_cast<real_long>( loc[0][n] );
        y0 = static_cast<real_long>( loc[1][n] );
        z0 = static_cast<real_long>( loc[2][n] );

        // rel prtcl index; assuming dx = 1; tile coordinates
        iff = D >= 1 ? static_cast<int>(floor( x0 - mins[0] ) ) : 0;
        jff = D >= 2 ? static_cast<int>(floor( y0 - mins[1] ) ) : 0;
        kff = D >= 3 ? static_cast<int>(floor( z0 - mins[2] ) ) : 0;

        // limit to 0 Nx just in case
        if(D >= 1) iff = limit(iff, 0, tile.mesh_lengths[0]-1);
        if(D >= 2) jff = limit(jff, 0, tile.mesh_lengths[1]-1);
        if(D >= 3) kff = limit(kff, 0, tile.mesh_lengths[2]-1);

        // update rho arrays
        yee.rho(iff,jff,kff) += mass;


        // full prtcl index; assuming dx = 1; global grid coordinates
  	    i = D >= 1 ? static_cast<int>(floor( x0 ) ) : 0;
  	    j = D >= 2 ? static_cast<int>(floor( y0 ) ) : 0;
  	    k = D >= 3 ? static_cast<int>(floor( z0 ) ) : 0;

        // reduce by a factor of stride; floating point arithmetics takes care of rounding
        if(D >= 1) i = limit(i/stride, 0, nx-1);
        if(D >= 2) j = limit(j/stride, 0, ny-1);
        if(D >= 3) k = limit(k/stride, 0, nz-1);

        u0 = static_cast<real_long>(vel[0][n]);
        v0 = static_cast<real_long>(vel[1][n]);
        w0 = static_cast<real_long>(vel[2][n]);

        gam = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);

        //--------------------------------------------------
        // next, physics quantities

        // number density 
        if(ispc == 0) dense(i,j,k) += 1.; 
        if(ispc == 1) densp(i,j,k) += 1.; 

        // bulk flows
        if(ispc == 0) {
          Vxe(i,j,k) += u0/gam; 
          Vye(i,j,k) += v0/gam;
          Vze(i,j,k) += w0/gam;
        }
        if(ispc == 1) {
          Vxp(i,j,k) += u0/gam; 
          Vyp(i,j,k) += v0/gam;
          Vzp(i,j,k) += w0/gam;
        }

        // momentum density (flux of mass)
        // momx(i,j,k)    += u0*mass;
        // momy(i,j,k)    += v0*mass;
        // momz(i,j,k)    += w0*mass;
          
        // pressure (flux of momentum)
        pressx(i,j,k)  += u0*u0*mass/gam;
        pressy(i,j,k)  += v0*v0*mass/gam;
        pressz(i,j,k)  += w0*w0*mass/gam;

        // off-diagonal shear terms 
        shearxy(i,j,k) += u0*v0*mass/gam;
        shearxz(i,j,k) += u0*w0*mass/gam;
        shearyz(i,j,k) += v0*w0*mass/gam;

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
template class h5io::PicMomentsWriter<2>;
template class h5io::PicMomentsWriter<3>;
