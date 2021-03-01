#include "field_slices.h"
#include <mpi4cpp/mpi.h>

#include "../../tools/ezh5/src/ezh5.hpp"
#include "../../em-fields/tile.h"


using namespace mpi4cpp;
using ezh5::File;


void h5io::FieldSliceWriter::read_tiles(
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

    // tile limits taking into account 0 collapsing dimensions
    int nxt, nyt, nzt;
    nxt = (int)yee.Nx/stride;
    nyt = (int)yee.Ny/stride;
    nzt = (int)yee.Nz/stride;

    nxt = nxt == 0 ? 1 : nxt;
    nyt = nyt == 0 ? 1 : nyt;
    nzt = nzt == 0 ? 1 : nzt;

    // starting location
    int i0 = nxt*std::get<0>(index);
    int j0 = nyt*std::get<1>(index);
    int k0 = nzt*std::get<2>(index);


    // x-y plane, z = 0
    if(mode == 0) {
      if(k0 != 0) continue;

      // field quantities 
      for(int js=0; js<nyt; js++) 
      for(int is=0; is<nxt; is++) {
        ex(i0+is, j0+js, 0) = yee.ex( is*stride, js*stride, 0);
        ey(i0+is, j0+js, 0) = yee.ey( is*stride, js*stride, 0);
        ez(i0+is, j0+js, 0) = yee.ez( is*stride, js*stride, 0);

        bx(i0+is, j0+js, 0) = yee.bx( is*stride, js*stride, 0);
        by(i0+is, j0+js, 0) = yee.by( is*stride, js*stride, 0);
        bz(i0+is, j0+js, 0) = yee.bz( is*stride, js*stride, 0);
      }

      // densities
      for(int js=0; js<nyt; js++) 
      for(int jstride=0; jstride < stride; jstride++) 
      for(int is=0; is<nxt; is++) 
      for(int istride=0; istride < stride; istride++) {
        jx(i0+is, j0+js, 0) += yee.jx( is*stride+istride, js*stride+jstride, 0);
        jy(i0+is, j0+js, 0) += yee.jy( is*stride+istride, js*stride+jstride, 0);
        jz(i0+is, j0+js, 0) += yee.jz( is*stride+istride, js*stride+jstride, 0);
        rh(i0+is, j0+js, 0) += yee.rho(is*stride+istride, js*stride+jstride, 0);
      }


    // x-z plane, y = 0
    } else if(mode == 1) {
      if(j0 != 0) continue;

      // field quantities
      for(int ks=0; ks<nzt; ks++) 
      for(int is=0; is<nxt; is++) {
        ex(i0+is, k0+ks, 0) = yee.ex( is*stride, 0, ks*stride);
        ey(i0+is, k0+ks, 0) = yee.ey( is*stride, 0, ks*stride);
        ez(i0+is, k0+ks, 0) = yee.ez( is*stride, 0, ks*stride);

        bx(i0+is, k0+ks, 0) = yee.bx( is*stride, 0, ks*stride);
        by(i0+is, k0+ks, 0) = yee.by( is*stride, 0, ks*stride);
        bz(i0+is, k0+ks, 0) = yee.bz( is*stride, 0, ks*stride);
      }

      // densities
      for(int ks=0; ks<nzt; ks++) 
      for(int kstride=0; kstride < stride; kstride++) 
      for(int is=0; is<nxt; is++) 
      for(int istride=0; istride < stride; istride++) {
        jx(i0+is, k0+ks, 0) += yee.jx( is*stride+istride, 0, ks*stride+kstride);
        jy(i0+is, k0+ks, 0) += yee.jy( is*stride+istride, 0, ks*stride+kstride);
        jz(i0+is, k0+ks, 0) += yee.jz( is*stride+istride, 0, ks*stride+kstride);
        rh(i0+is, k0+ks, 0) += yee.rho(is*stride+istride, 0, ks*stride+kstride);
      }

    // y-z plane, x = 0
    } else if(mode == 2) {
      if(i0 != 0) continue;

      // field quantities; just downsample by hopping with stride
      for(int ks=0; ks<nzt; ks++) 
      for(int js=0; js<nyt; js++) {
        ex(j0+js, k0+ks, 0) = yee.ex(0, js*stride, ks*stride);
        ey(j0+js, k0+ks, 0) = yee.ey(0, js*stride, ks*stride);
        ez(j0+js, k0+ks, 0) = yee.ez(0, js*stride, ks*stride);

        bx(j0+js, k0+ks, 0) = yee.bx(0, js*stride, ks*stride);
        by(j0+js, k0+ks, 0) = yee.by(0, js*stride, ks*stride);
        bz(j0+js, k0+ks, 0) = yee.bz(0, js*stride, ks*stride);
      }

      // densities; these quantities we average over the volume
      for(int ks=0; ks<nzt; ks++) 
      for(int kstride=0; kstride < stride; kstride++) 
      for(int js=0; js<nyt; js++) 
      for(int jstride=0; jstride < stride; jstride++) {
        jx(j0+js, k0+ks, 0) += yee.jx( 0, js*stride+jstride, ks*stride+kstride);
        jy(j0+js, k0+ks, 0) += yee.jy( 0, js*stride+jstride, ks*stride+kstride);
        jz(j0+js, k0+ks, 0) += yee.jz( 0, js*stride+jstride, ks*stride+kstride);
        rh(j0+js, k0+ks, 0) += yee.rho(0, js*stride+jstride, ks*stride+kstride);
      }
    }
  } // tiles
}

bool h5io::FieldSliceWriter::write(
    corgi::Grid<3>& grid, int lap)
{
  read_tiles(grid);
  mpi_reduce_snapshots(grid);

  if( grid.comm.rank() == 0 ) {

    // build filename

    std::string mode_name;
    if(mode == 0) {
      mode_name = "xy";
    } else if(mode == 1) {
      mode_name = "xz";
    } else if(mode == 2) {
      mode_name = "yz";
    }

    std::string full_filename =
      fname +
      +"/"+
      file_name +
      "-" +
      mode_name +
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


