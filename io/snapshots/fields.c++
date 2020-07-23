#include "fields.h"
#include <mpi4cpp/mpi.h>

#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../em-fields/tile.h"


using namespace mpi4cpp;
using ezh5::File;


template<>
inline void h5io::FieldsWriter<2>::read_tiles(
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
inline void h5io::FieldsWriter<3>::read_tiles(
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
inline bool h5io::FieldsWriter<D>::write(
    corgi::Grid<D>& grid, int lap)
{
  read_tiles(grid);
  mpi_reduce_snapshots(grid);

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
//template class h5io::FieldsWriter<1>;
template class h5io::FieldsWriter<2>;
template class h5io::FieldsWriter<3>;
