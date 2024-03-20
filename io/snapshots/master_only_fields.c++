#include "master_only_fields.h"
#include <mpi4cpp/mpi.h>

#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../emf/tile.h"


using namespace mpi4cpp;
using ezh5::File;


template<>
inline void h5io::MasterFieldsWriter<3>::read_tiles(
    corgi::Grid<3>& grid)
{
  // this function is never used with this class
  assert(false);
}


template<>
inline void h5io::MasterFieldsWriter<3>::read_tile_feature(
    corgi::Grid<3>& grid,
    uint64_t cid,
    int ifea
    )
{
  auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile( cid ));
  auto& yee = tile.get_yee();
    
  // clear buffer before additive variables
  sbuf[0].clear();

  if(ifea <= 5) {

    // field quantities; just downsample by hopping with stride
    for(int ks=0; ks<nzM; ks++) 
    for(int js=0; js<nyM; js++) 
    for(int is=0; is<nxM; is++) {
      if(ifea == 0) sbuf[0](is, js, ks) = yee.ex( is*stride, js*stride, ks*stride);
      if(ifea == 1) sbuf[0](is, js, ks) = yee.ey( is*stride, js*stride, ks*stride);
      if(ifea == 2) sbuf[0](is, js, ks) = yee.ez( is*stride, js*stride, ks*stride);

      if(ifea == 3) sbuf[0](is, js, ks) = yee.bx( is*stride, js*stride, ks*stride);
      if(ifea == 4) sbuf[0](is, js, ks) = yee.by( is*stride, js*stride, ks*stride);
      if(ifea == 5) sbuf[0](is, js, ks) = yee.bz( is*stride, js*stride, ks*stride);
    }

  } else {

    // densities; these quantities we average over the volume
    for(int ks=0; ks<nzM; ks++) 
    for(int kstride=0; kstride < stride; kstride++) 
    for(int js=0; js<nyM; js++) 
    for(int jstride=0; jstride < stride; jstride++) 
    for(int is=0; is<nxM; is++) 
    for(int istride=0; istride < stride; istride++) {
      if(ifea == 6) sbuf[0](is, js, ks) += yee.jx( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      if(ifea == 7) sbuf[0](is, js, ks) += yee.jy( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      if(ifea == 8) sbuf[0](is, js, ks) += yee.jz( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      if(ifea == 9) sbuf[0](is, js, ks) += yee.rho(is*stride+istride, js*stride+jstride, ks*stride+kstride);
    }

  }

  return;
}



template<>
inline void h5io::MasterFieldsWriter<3>::mpi_reduce_snapshots(
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
  int n_tiles = nx_tile*ny_tile*nz_tile;

  // sync everyone before going into the loop
  comm.barrier();

  for(uint64_t cid=0; cid<n_tiles; cid++) {

    int msg_rank = -1;
    if( grid.is_local(cid) ) msg_rank = rank;
    bool my_msg  = msg_rank == rank; // check if message rank matches the local one 

    for(int ifea=0; ifea<10; ifea++){

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
inline bool h5io::MasterFieldsWriter<3>::write(
    corgi::Grid<3>& grid, int lap)
{

  //--------------------------------------------------
  // allocate full array for master
  if( grid.comm.rank() == 0) {

    // allocate if first time 
    if(!master_is_initialized){
      for(size_t i=0; i<10; i++) arrs.emplace_back(nx, ny, nz);
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
    
    file["ex"] = arrs[0].serialize();
    file["ey"] = arrs[1].serialize();
    file["ez"] = arrs[2].serialize();

    file["bx"] = arrs[3].serialize();
    file["by"] = arrs[4].serialize();
    file["bz"] = arrs[5].serialize();

    file["jx"] = arrs[6].serialize();
    file["jy"] = arrs[7].serialize();
    file["jz"] = arrs[8].serialize();

    file["rho"]= arrs[9].serialize();
    
  }

  return true;
}


//--------------------------------------------------
// explicit template class instantiations
//template class h5io::MasterFieldsWriter<1>;
//template class h5io::MasterFieldsWriter<2>;
template class h5io::MasterFieldsWriter<3>;
