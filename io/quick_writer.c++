#include "quick_writer.h"

#include <mpi4cpp/mpi.h>
#include "../tools/mesh.h"
#include "../em-fields/tile.h"
#include "../tools/ezh5/src/ezh5.hpp"
#include "../tools/fastlog.h"
#include "namer.h"

using namespace mpi4cpp;
using ezh5::File;


template<>
inline void h5io::QuickWriter<2>::read_tiles(
    corgi::Grid<2>& grid)
{

  // clear target arrays
  for(auto& arr : arrs) arr.clear();

  // target arrays
  auto& ex = arrs[0];
  auto& ey = arrs[1];
  auto& ez = arrs[2];
  auto& bx = arrs[3];
  auto& by = arrs[4];
  auto& bz = arrs[5];
  auto& jx = arrs[6];
  auto& jy = arrs[7];
  auto& jz = arrs[8];
  auto& rh = arrs[9];


  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<fields::Tile<2>&>(grid.get_tile( cid ));
    auto& yee = tile.get_yee();

    // get arrays
    auto index = expand_indices( &tile );

    // starting location
    int i0 = (yee.Nx/stride)*std::get<0>(index);
    int j0 = (yee.Ny/stride)*std::get<1>(index);
    int k0 = 0; //(yee.Nz/stride)*std::get<2>(index);

    // tile limits taking into account 0 collapsing dimensions
    int nxt, nyt;
    nxt = (int)yee.Nx/stride;
    nyt = (int)yee.Ny/stride;

    nxt = nxt == 0 ? 1 : nxt;
    nyt = nyt == 0 ? 1 : nyt;

    // copy tile patch by stride hopping
    int ks = 0;
    int kstride = 0;

    for(int js=0; js<nyt; js++) {
      for(int jstride=0; jstride < stride; jstride++) {
        for(int is=0; is<nxt; is++) {

          // field quantities; no integration
          if (jstride==0) {
            ex(i0+is, j0+js, k0+ks) += yee.ex( is*stride, js*stride, ks*stride);
            ey(i0+is, j0+js, k0+ks) += yee.ey( is*stride, js*stride, ks*stride);
            ez(i0+is, j0+js, k0+ks) += yee.ez( is*stride, js*stride, ks*stride);

            bx(i0+is, j0+js, k0+ks) += yee.bx( is*stride, js*stride, ks*stride);
            by(i0+is, j0+js, k0+ks) += yee.by( is*stride, js*stride, ks*stride);
            bz(i0+is, j0+js, k0+ks) += yee.bz( is*stride, js*stride, ks*stride);
          }

          // densities; these quantities we need to integrate over stride
          for(int istride=0; istride < stride; istride++) {
            jx(i0+is, j0+js, k0+ks) += yee.jx( is*stride+istride, js*stride+jstride, ks*stride+kstride);
            jy(i0+is, j0+js, k0+ks) += yee.jy( is*stride+istride, js*stride+jstride, ks*stride+kstride);
            jz(i0+is, j0+js, k0+ks) += yee.jz( is*stride+istride, js*stride+jstride, ks*stride+kstride);
            rh(i0+is, j0+js, k0+ks) += yee.rho(is*stride+istride, js*stride+jstride, ks*stride+kstride);
          }

        }
      }
    }

  } // tiles
}


template<>
inline void h5io::QuickWriter<3>::read_tiles(
    corgi::Grid<3>& grid)
{

  // clear target arrays
  for(auto& arr : arrs) arr.clear();

  // target arrays
  auto& ex = arrs[0];
  auto& ey = arrs[1];
  auto& ez = arrs[2];
  auto& bx = arrs[3];
  auto& by = arrs[4];
  auto& bz = arrs[5];
  auto& jx = arrs[6];
  auto& jy = arrs[7];
  auto& jz = arrs[8];
  auto& rh = arrs[9];


  // read my local tiles
  for(auto cid : grid.get_local_tiles() ){
    auto& tile = dynamic_cast<fields::Tile<3>&>(grid.get_tile( cid ));
    auto& yee = tile.get_yee();

    // get arrays
    auto index = expand_indices( &tile );

    // starting location
    int i0 = (yee.Nx/stride)*std::get<0>(index);
    int j0 = (yee.Ny/stride)*std::get<1>(index);
    int k0 = (yee.Nz/stride)*std::get<2>(index);

    // tile limits taking into account 0 collapsing dimensions
    int nxt, nyt, nzt;
    nxt = (int)yee.Nx/stride;
    nyt = (int)yee.Ny/stride;
    nzt = (int)yee.Nz/stride;

    nxt = nxt == 0 ? 1 : nxt;
    nyt = nyt == 0 ? 1 : nyt;
    nzt = nzt == 0 ? 1 : nzt;

    // copy tile patch by stride hopping; either downsample or average

    // field quantities; just downsample by hopping with stride
    for(int ks=0; ks<nzt; ks++) 
    for(int js=0; js<nyt; js++) 
    for(int is=0; is<nxt; is++) {
      ex(i0+is, j0+js, k0+ks) = yee.ex( is*stride, js*stride, ks*stride);
      ey(i0+is, j0+js, k0+ks) = yee.ey( is*stride, js*stride, ks*stride);
      ez(i0+is, j0+js, k0+ks) = yee.ez( is*stride, js*stride, ks*stride);

      bx(i0+is, j0+js, k0+ks) = yee.bx( is*stride, js*stride, ks*stride);
      by(i0+is, j0+js, k0+ks) = yee.by( is*stride, js*stride, ks*stride);
      bz(i0+is, j0+js, k0+ks) = yee.bz( is*stride, js*stride, ks*stride);
    }

    // densities; these quantities we average over the volume
    for(int ks=0; ks<nzt; ks++) 
    for(int kstride=0; kstride < stride; kstride++) 
    for(int js=0; js<nyt; js++) 
    for(int jstride=0; jstride < stride; jstride++) 
    for(int is=0; is<nxt; is++) 
    for(int istride=0; istride < stride; istride++) {
      jx(i0+is, j0+js, k0+ks) += yee.jx( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      jy(i0+is, j0+js, k0+ks) += yee.jy( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      jz(i0+is, j0+js, k0+ks) += yee.jz( is*stride+istride, js*stride+jstride, ks*stride+kstride);
      rh(i0+is, j0+js, k0+ks) += yee.rho(is*stride+istride, js*stride+jstride, ks*stride+kstride);
    }

  } // tiles
}







template<size_t D>
inline void h5io::QuickWriter<D>::mpi_reduce_snapshots(
    corgi::Grid<D>& grid)
{
  /* based on https://gist.github.com/rmcgibbo/7178576
   */

  for(auto& arr : rbuf) arr.clear();

  int tag = 0;
  const int size = grid.comm.size(); //MPI::COMM_WORLD.Get_size();
  const int rank = grid.comm.rank(); //MPI::COMM_WORLD.Get_rank();
  const int lastpower = 1 << fastlog2(size);

  // each of the ranks greater than the last power of 2 less than size
  // need to downshift their data, since the binary tree reduction below
  // only works when N is a power of two.
  
  std::vector<mpi::request> reqs;
  for (int i = lastpower; i < size; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        reqs.push_back( grid.comm.isend(i-lastpower, tag+els, arrs[els].data(), arrs[els].size()) );
      }
    }
  }

  for (int i = 0; i < size-lastpower; i++) {
    if (rank == i) {
      for(size_t els=0; els<arrs.size(); els++) {
        grid.comm.recv(i+lastpower, tag+els, rbuf[els].data(), rbuf[els].size());
        arrs[els] += rbuf[els];
      }
    }
  }
  mpi::wait_all(reqs.begin(), reqs.end());


  for (int d = 0; d < fastlog2(lastpower); d++) {
    for (int k = 0; k < lastpower; k += 1 << (d + 1)) {
      const int receiver = k;
      const int sender = k + (1 << d);
      if (rank == receiver) {

        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.recv(sender, tag+els, rbuf[els].data(), rbuf[els].size());
          arrs[els] += rbuf[els];
        }
      }
      else if (rank == sender) {
        for(size_t els=0; els<arrs.size(); els++) {
          grid.comm.send(receiver, tag+els, arrs[els].data(), arrs[els].size());
        }
      }
    }
  }


}


template<size_t D>
inline bool h5io::QuickWriter<D>::write(
    corgi::Grid<D>& grid, int lap)
{
  read_tiles(grid);
  mpi_reduce_snapshots(grid);

  if( grid.comm.rank() == 0 ) {
    // build filename
    std::string full_filename = 
      fname + 
      "/flds" +
      /*std::to_string(grid.comm.rank()) + */
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
    file["ex"] = arrs[0].mat;
    file["ey"] = arrs[1].mat;
    file["ez"] = arrs[2].mat;

    file["bx"] = arrs[3].mat;
    file["by"] = arrs[4].mat;
    file["bz"] = arrs[5].mat;

    file["jx"] = arrs[6].mat;
    file["jy"] = arrs[7].mat;
    file["jz"] = arrs[8].mat;

    file["rho"]= arrs[9].mat;
  }

  return true;
}






//--------------------------------------------------
// explicit template class instantiations
//template class h5io::QuickWriter<1>;
template class h5io::QuickWriter<2>;
template class h5io::QuickWriter<3>;
