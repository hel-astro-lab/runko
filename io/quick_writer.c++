#include "quick_writer.h"

#include "../tools/mesh.h"
#include "../em-fields/tile.h"
#include "../tools/ezh5/src/ezh5.hpp"
#include "namer.h"


using ezh5::File;
template<size_t D>
inline void h5io::QuickWriter<D>::read_tiles(
    corgi::Node<D>& grid)
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
  for(auto cid : grid.get_tile_ids() ){
    auto& tile = dynamic_cast<fields::Tile<D>&>(grid.get_tile( cid ));
    auto& yee = tile.get_yee();

    // get arrays
    // strides: 
    //  average or just closest value
    auto index = expand_indices( &tile );

    // starting location
    int i0 = (yee.Nx/stride)*std::get<0>(index);
    int j0 = (yee.Ny/stride)*std::get<1>(index);
    int k0 = (yee.Nz/stride)*std::get<2>(index);

    // copy tile patch by stride hopping
    for(int ks=0; ks<(int)yee.Nz/stride; ks++ ) {
      for(int kstride=0; kstride < stride; kstride++) {

        for(int js=0; js<(int)yee.Ny/stride; js++) {
          for(int jstride=0; jstride < stride; jstride++) {

            for(int is=0; is<(int)yee.Nx/stride; is++) {
              for(int istride=0; istride < stride; istride++) {
                ex(i0+is, j0+js, k0+ks) += yee.ex( is+istride, js+jstride, ks+stride);
                ey(i0+is, j0+js, k0+ks) += yee.ey( is+istride, js+jstride, ks+stride);
                ez(i0+is, j0+js, k0+ks) += yee.ez( is+istride, js+jstride, ks+stride);
                                                                                     
                bx(i0+is, j0+js, k0+ks) += yee.bx( is+istride, js+jstride, ks+stride);
                by(i0+is, j0+js, k0+ks) += yee.by( is+istride, js+jstride, ks+stride);
                bz(i0+is, j0+js, k0+ks) += yee.bz( is+istride, js+jstride, ks+stride);
                                                                                     
                jx(i0+is, j0+js, k0+ks) += yee.jx( is+istride, js+jstride, ks+stride);
                jy(i0+is, j0+js, k0+ks) += yee.jy( is+istride, js+jstride, ks+stride);
                jz(i0+is, j0+js, k0+ks) += yee.jz( is+istride, js+jstride, ks+stride);
                                                                                     
                rh(i0+is, j0+js, k0+ks) += yee.rho(is+istride, js+jstride, ks+stride);
              }
            }
          }
        }
      }
    }
  }


}



template<size_t D>
inline bool h5io::QuickWriter<D>::write(
    corgi::Node<D>& grid, int lap)
{
  read_tiles(grid);





  // build filename
  std::string full_filename = 
    fname + 
    "-" +
    std::to_string(grid.comm.rank()) + 
    "_" +
    std::to_string(lap) +
    extension;


  // open file and write
  file(full_filename, H5F_ACC_TRUNC);

  file["ex"] = arr[0].serialize();
  file["ey"] = arr[1].serialize();
  file["ez"] = arr[2].serialize();

  file["bx"] = arr[3].serialize();
  file["by"] = arr[4].serialize();
  file["bz"] = arr[5].serialize();

  file["jx"] = arr[6].serialize();
  file["jy"] = arr[7].serialize();
  file["jz"] = arr[8].serialize();

  file["rho"]= arr[9].serialize();

  return true;
}






//--------------------------------------------------
// explicit template class instantiations
template class h5io::QuickWriter<1>;
template class h5io::QuickWriter<2>;
