#include <mpi4cpp/mpi.h>

#include "io/snapshots/field_slices.h"
#include "external/ezh5/src/ezh5.hpp"
#include "core/emf/tile.h"


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
    auto& tile = dynamic_cast<emf::Tile<3>&>(grid.get_tile( cid ));
    auto& gs = tile.get_grids();

    // get arrays
    auto index = expand_indices( &tile );

    // tile limits taking into account 0 collapsing dimensions
    int nxt, nyt, nzt;
    nxt = (int)gs.Nx/stride;
    nyt = (int)gs.Ny/stride;
    nzt = (int)gs.Nz/stride;

    nxt = nxt == 0 ? 1 : nxt;
    nyt = nyt == 0 ? 1 : nyt;
    nzt = nzt == 0 ? 1 : nzt;

    // starting location
    int i0 = nxt*std::get<0>(index);
    int j0 = nyt*std::get<1>(index);
    int k0 = nzt*std::get<2>(index);

    bool ind_inside_x = tile.mins[0] <= ind && ind < tile.maxs[0];
    bool ind_inside_y = tile.mins[1] <= ind && ind < tile.maxs[1];
    bool ind_inside_z = tile.mins[2] <= ind && ind < tile.maxs[2];

    // x-y plane, z = ind
    if(mode == 0) {

      //if(k0 != 0) continue;
      if(!ind_inside_z) continue;

      int I = ind - tile.mins[2]; // location on the tile

      // field quantities 
      for(int js=0; js<nyt; js++) 
      for(int is=0; is<nxt; is++) {
        ex(i0+is, j0+js, 0) = gs.ex( is*stride, js*stride, I);
        ey(i0+is, j0+js, 0) = gs.ey( is*stride, js*stride, I);
        ez(i0+is, j0+js, 0) = gs.ez( is*stride, js*stride, I);

        bx(i0+is, j0+js, 0) = gs.bx( is*stride, js*stride, I);
        by(i0+is, j0+js, 0) = gs.by( is*stride, js*stride, I);
        bz(i0+is, j0+js, 0) = gs.bz( is*stride, js*stride, I);
      }

      // densities
      for(int js=0; js<nyt; js++) 
      for(int jstride=0; jstride < stride; jstride++) 
      for(int is=0; is<nxt; is++) 
      for(int istride=0; istride < stride; istride++) {
        jx(i0+is, j0+js, 0) += gs.jx( is*stride+istride, js*stride+jstride, I);
        jy(i0+is, j0+js, 0) += gs.jy( is*stride+istride, js*stride+jstride, I);
        jz(i0+is, j0+js, 0) += gs.jz( is*stride+istride, js*stride+jstride, I);
        rh(i0+is, j0+js, 0) += gs.rho(is*stride+istride, js*stride+jstride, I);
      }


    // x-z plane, y = 0
    } else if(mode == 1) {

      //if(j0 != 0) continue;
      if(!ind_inside_y) continue;

      int I = ind - tile.mins[1]; // location on the tile

      // field quantities
      for(int ks=0; ks<nzt; ks++) 
      for(int is=0; is<nxt; is++) {
        ex(i0+is, k0+ks, 0) = gs.ex( is*stride, I, ks*stride);
        ey(i0+is, k0+ks, 0) = gs.ey( is*stride, I, ks*stride);
        ez(i0+is, k0+ks, 0) = gs.ez( is*stride, I, ks*stride);

        bx(i0+is, k0+ks, 0) = gs.bx( is*stride, I, ks*stride);
        by(i0+is, k0+ks, 0) = gs.by( is*stride, I, ks*stride);
        bz(i0+is, k0+ks, 0) = gs.bz( is*stride, I, ks*stride);
      }

      // densities
      for(int ks=0; ks<nzt; ks++) 
      for(int kstride=0; kstride < stride; kstride++) 
      for(int is=0; is<nxt; is++) 
      for(int istride=0; istride < stride; istride++) {
        jx(i0+is, k0+ks, 0) += gs.jx( is*stride+istride, I, ks*stride+kstride);
        jy(i0+is, k0+ks, 0) += gs.jy( is*stride+istride, I, ks*stride+kstride);
        jz(i0+is, k0+ks, 0) += gs.jz( is*stride+istride, I, ks*stride+kstride);
        rh(i0+is, k0+ks, 0) += gs.rho(is*stride+istride, I, ks*stride+kstride);
      }

    // y-z plane, x = 0
    } else if(mode == 2) {

      //if(i0 != 0) continue;
      if(!ind_inside_x) continue;

      int I = ind - tile.mins[0]; // location on the tile

      // field quantities; just downsample by hopping with stride
      for(int ks=0; ks<nzt; ks++) 
      for(int js=0; js<nyt; js++) {
        ex(j0+js, k0+ks, 0) = gs.ex(I, js*stride, ks*stride);
        ey(j0+js, k0+ks, 0) = gs.ey(I, js*stride, ks*stride);
        ez(j0+js, k0+ks, 0) = gs.ez(I, js*stride, ks*stride);

        bx(j0+js, k0+ks, 0) = gs.bx(I, js*stride, ks*stride);
        by(j0+js, k0+ks, 0) = gs.by(I, js*stride, ks*stride);
        bz(j0+js, k0+ks, 0) = gs.bz(I, js*stride, ks*stride);
      }

      // densities; these quantities we average over the volume
      for(int ks=0; ks<nzt; ks++) 
      for(int kstride=0; kstride < stride; kstride++) 
      for(int js=0; js<nyt; js++) 
      for(int jstride=0; jstride < stride; jstride++) {
        jx(j0+js, k0+ks, 0) += gs.jx( I, js*stride+jstride, ks*stride+kstride);
        jy(j0+js, k0+ks, 0) += gs.jy( I, js*stride+jstride, ks*stride+kstride);
        jz(j0+js, k0+ks, 0) += gs.jz( I, js*stride+jstride, ks*stride+kstride);
        rh(j0+js, k0+ks, 0) += gs.rho(I, js*stride+jstride, ks*stride+kstride);
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


