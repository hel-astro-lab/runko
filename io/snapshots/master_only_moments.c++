#include "master_only_moments.h"
#include <mpi4cpp/mpi.h>

#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../pic/particle.h"
#include "../../pic/tile.h"
#include "../../tools/signum.h"
#include "../../tools/limit.h"


using namespace mpi4cpp;
using ezh5::File;


template<>
inline void h5io::MasterPicMomentsWriter<3>::read_tiles(
    corgi::Grid<3>& grid)
{
  // this function is never used with this class
  assert(false);
}


template<size_t D>
inline void h5io::MasterPicMomentsWriter<D>::read_tile_feature(
    corgi::Grid<D>& grid,
    uint64_t cid,
    int ifea
    )
{
  auto& tile = dynamic_cast<pic::Tile<D>&>(grid.get_tile( cid ));
  auto& yee = tile.get_grids();
    
  auto mins = tile.mins;
  auto maxs = tile.maxs;


  // clear buffer before additive variables
  sbuf[0].clear();
  if(ifea==0) yee.rho.clear();


  // local variables
  double gam, xene, mass, wgt;
  double x0, y0, z0;
  double u0, v0, w0;
  int nparts;
  int i,j,k;
  int iff,jff,kff;



  // loop over species
  for (int ispc=0; ispc<tile.Nspecies(); ispc++) {
    auto& container = tile.get_container(ispc);
    mass = container.m; // species mass
    nparts = container.size();

    if(nparts <= 0) continue; // skip zero containers

    float_p* loc[3];
    for( i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

    float_p* vel[3];
    for( i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

    float_p* ch;
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

      gam  = sqrt(1.0 + u0*u0 + v0*v0 + w0*w0);
      xene = sqrt(      u0*u0 + v0*v0 + w0*w0);


      //--------------------------------------------------
      if(ifea == 0) { // only first feature checks these tests and adds rho

        // capture NaNs
        bool isnan1 = std::isnan(x0);
        bool isnan2 = std::isnan(y0);
        bool isnan3 = std::isnan(z0);
        bool isnan4 = std::isnan(wgt);

        // check that particles are communicated properly
        bool flagx = D>= 1 ? (x0-mins[0] >= -3.0) && (x0 <= maxs[0] +2.0 ) : true;
        bool flagy = D>= 2 ? (y0-mins[1] >= -3.0) && (y0 <= maxs[1] +2.0 ) : true;
        bool flagz = D>= 3 ? (z0-mins[2] >= -3.0) && (z0 <= maxs[2] +2.0 ) : true;

        if( !flagx || !flagy || !flagz || 
            isnan1 || isnan2 || isnan3 || isnan4) {
          std::cerr << "ERR IN MOM:" << std::endl;
          std::cerr << " ispc: " << ispc;
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
          //assert(false);
        }

        // rel prtcl index; assuming dx = 1; tile coordinates
        // limit to 0 Nx-1 just in case to avoid crashes
        iff = D >= 1 ? limit( floor(x0-mins[0]), -3., maxs[0]-mins[0] +2.) : 0;
        jff = D >= 2 ? limit( floor(y0-mins[1]), -3., maxs[1]-mins[1] +2.) : 0;
        kff = D >= 3 ? limit( floor(z0-mins[2]), -3., maxs[2]-mins[2] +2.) : 0;

        // update rho arrays; this is interpreted as mass density
        yee.rho(iff,jff,kff) += mass*wgt;

      //--------------------------------------------------
      }


      // full prtcl index; assuming dx = 1; global grid coordinates reduce by a factor of stride 
      // NOTE in tile-origin based coordinate system; transformed back to global indexing (by root) 
      //      after mpi messages are done.
      i = D >= 1 ? limit( floor( (x0-mins[0])/stride), 0.0, double(nxM)-1.0) : 0;
      j = D >= 2 ? limit( floor( (y0-mins[1])/stride), 0.0, double(nyM)-1.0) : 0;
      k = D >= 3 ? limit( floor( (z0-mins[2])/stride), 0.0, double(nzM)-1.0) : 0;

      // number density 
      if(ifea == 0 && ispc == 0) sbuf[0](i,j,k) += wgt; 
      if(ifea == 1 && ispc == 1) sbuf[0](i,j,k) += wgt; 
      if(ifea == 2 && ispc == 2) sbuf[0](i,j,k) += wgt*xene;  // energy density for photons

      // bulk flows
      if(ifea == 3 && ispc == 0) sbuf[0](i,j,k) += wgt*u0/gam; // Vxe
      if(ifea == 4 && ispc == 0) sbuf[0](i,j,k) += wgt*v0/gam; // Vye
      if(ifea == 5 && ispc == 0) sbuf[0](i,j,k) += wgt*w0/gam; // Vze

      if(ifea == 6 && ispc == 1) sbuf[0](i,j,k) += wgt*u0/gam; // Vxp
      if(ifea == 7 && ispc == 1) sbuf[0](i,j,k) += wgt*v0/gam; // Vyp
      if(ifea == 8 && ispc == 1) sbuf[0](i,j,k) += wgt*w0/gam; // Vzp


      // pressure (flux of momentum)
      if(ifea == 9             ) sbuf[0](i,j,k) += wgt*u0*u0*mass/gam; // pressx
      if(ifea ==10             ) sbuf[0](i,j,k) += wgt*v0*v0*mass/gam; // pressy
      if(ifea ==11             ) sbuf[0](i,j,k) += wgt*w0*w0*mass/gam; // pressz
        
      // off-diagonal shear terms 
      if(ifea ==12             ) sbuf[0](i,j,k) += wgt*u0*v0*mass/gam; // shearxy
      if(ifea ==13             ) sbuf[0](i,j,k) += wgt*v0*w0*mass/gam; // shearxz
      if(ifea ==14             ) sbuf[0](i,j,k) += wgt*v0*w0*mass/gam; // shearyz

    }
  }

  // FIXME cant normalize Ve's and Vp's because we dont store dense/densp
  //if( ifea == 3 ||
  //    ifea == 4 ||
  //    ifea == 5 ||
  //    //
  //    ifea == 6 ||
  //    ifea == 7 ||
  //    ifea == 8 ){

  //  for(int k=0; k<nzM; k++)
  //  for(int j=0; j<nyM; j++)
  //  for(int i=0; i<nxM; i++) {

  //  }
  //}


  return;
}



template<>
inline void h5io::MasterPicMomentsWriter<3>::mpi_reduce_snapshots(
    corgi::Grid<3>& grid)
{
  // mpi communicator
  auto& comm = grid.comm;
     
  int rank = comm.rank(); // local rank
  bool is_master = rank == 0; 

  //int n_tiles = grid._mpi_grid.size();
  auto lens = grid.lens();
  int nx_tile = lens[0];
  int ny_tile = lens[1];
  int nz_tile = lens[2];
  size_t n_tiles = nx_tile*ny_tile*nz_tile;

  // sync everyone before going into the loop
  comm.barrier();

  for(size_t cid=0; cid<n_tiles; cid++) {

    int msg_rank = -1;
    if( grid.is_local(cid) ) msg_rank = rank;
    bool my_msg  = msg_rank == rank; // check if message rank matches the local one 

    for(int ifea=0; ifea<15; ifea++){

      //--------------------------------------------------
      // message sending

      // owner of the tile sends the message
      if(my_msg && !is_master) {

        // load feature into the sbuf
        read_tile_feature(grid, cid, ifea);

        // mpi send; using feature as the mpi tag to distinguish between packets
        comm.send(0, ifea, sbuf[0].data(),  nxM*nyM*nzM);

      } else if(my_msg && is_master) {
        // NOTE special branch for root. It does not send or receive
          
        // load feature into the sbuf
        read_tile_feature(grid, cid, ifea);

        // copy from sbuf to rbuf manually
        for(int ks=0; ks<nzM; ks++) 
        for(int js=0; js<nyM; js++) 
        for(int is=0; is<nxM; is++) {
          rbuf[0](is,js,ks) = sbuf[0](is,js,ks);
        }

      }

      //--------------------------------------------------
      // message receiving

      // master unpacks the message
      if(is_master) {
          
        // starting location
        auto index = grid.id2index(cid, grid.lens() );
        int i0 = nxM*std::get<0>(index);
        int j0 = nyM*std::get<1>(index);
        int k0 = nzM*std::get<2>(index);

        // mpi receive; using feature as the mpi tag to distinguish between packets
        // NOTE root does not send nor receive
        if(!my_msg && is_master) {
          comm.recv(msg_rank, ifea, rbuf[0].data(),  nxM*nyM*nzM);
        }

        // unpack to global grid
        for(int ks=0; ks<nzM; ks++) 
        for(int js=0; js<nyM; js++) 
        for(int is=0; is<nxM; is++) {
          arrs[ifea](i0+is, j0+js, k0+ks) = rbuf[0]( is, js, ks);
        }
      }

      // everyone waits, tile by tile. Slow but gives minimum memory footprint
      comm.barrier();
    }

  }

  return;
}




//template<size_t D>
template<>
inline bool h5io::MasterPicMomentsWriter<3>::write(
    corgi::Grid<3>& grid, int lap)
{

  //--------------------------------------------------
  // allocate full array for master
  if( grid.comm.rank() == 0) {

    // allocate if first time 
    if(!master_is_initialized){
      for(size_t i=0; i<15; i++) arrs.emplace_back(nx, ny, nz);
    }

    // otherwise clear target arrays
    for(auto& arr : arrs) arr.clear();
  }

  //--------------------------------------------------
  mpi_reduce_snapshots(grid);

  // root writes
  if( grid.comm.rank() == 0 ) {

    // build filename
    std::string full_filename =
      fname +
      +"/"+
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

    // avoid extra copy by using internal container reference;
    // this works because writer meshes don't have halos

    // NOTE index ordering is different than in multi-rank pic moment writer
    
    file["dense"] = arrs[0].serialize();
    file["densp"] = arrs[1].serialize();
    file["densx"] = arrs[2].serialize(); 

    file["Vxe"]   = arrs[3].serialize();
    file["Vye"]   = arrs[4].serialize();
    file["Vze"]   = arrs[5].serialize();

    file["Vxp"]   = arrs[6].serialize();
    file["Vyp"]   = arrs[7].serialize();
    file["Vzp"]   = arrs[8].serialize();

    file["pressx"]   = arrs[9].serialize();
    file["pressy"]   = arrs[10].serialize();
    file["pressz"]   = arrs[11].serialize();

    file["shearxy"]   = arrs[12].serialize();
    file["shearxz"]   = arrs[13].serialize();
    file["shearyz"]   = arrs[14].serialize();

  }

  return true;
}


//--------------------------------------------------
// explicit template class instantiations
template class h5io::MasterPicMomentsWriter<3>;
