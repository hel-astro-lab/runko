#include "quick_writer.h"

#include "../tools/mesh.h"
#include "../em-fields/tile.h"
#include "../tools/ezh5/src/ezh5.hpp"
#include "namer.h"


using ezh5::File;
//template<size_t D>
template<>
inline void h5io::QuickWriter<2>::read_tiles(
    corgi::Node<2>& grid)
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
    auto& tile = dynamic_cast<fields::Tile<2>&>(grid.get_tile( cid ));
    auto& yee = tile.get_yee();

    // get arrays
    // strides: 
    //  average or just closest value
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

    // copy tile patch by stride hopping
    //for(int ks=0; ks<nzt; ks++ ) {
    //  for(int kstride=0; kstride < stride; kstride++) {
    int ks = 0;
    int kstride = 0;
    k0 = 0;
        for(int js=0; js<nyt; js++) {
          for(int jstride=0; jstride < stride; jstride++) {

            for(int is=0; is<nxt; is++) {
              for(int istride=0; istride < stride; istride++) {
                //std::cout << "summing from "
                //  << " " << i0 << " " << is 
                //  << " " << j0 << " " << js 
                //  << " " << k0 << " " << ks 
                //  << " " << "to "
                //  << " " << is << " " << istride
                //  << " " << js << " " << jstride
                //  << " " << ks << " " << kstride
                //  << " val " << rh(i0+is, j0+js, k0+ks) << " += " << yee.rho(is+istride, js+jstride, ks+kstride)
                //  << std::endl;


                ex(i0+is, j0+js, k0+ks) += yee.ex( is+istride, js+jstride, ks+kstride);
                ey(i0+is, j0+js, k0+ks) += yee.ey( is+istride, js+jstride, ks+kstride);
                ez(i0+is, j0+js, k0+ks) += yee.ez( is+istride, js+jstride, ks+kstride);
                                                                                      
                bx(i0+is, j0+js, k0+ks) += yee.bx( is+istride, js+jstride, ks+kstride);
                by(i0+is, j0+js, k0+ks) += yee.by( is+istride, js+jstride, ks+kstride);
                bz(i0+is, j0+js, k0+ks) += yee.bz( is+istride, js+jstride, ks+kstride);
                                                                                      
                jx(i0+is, j0+js, k0+ks) += yee.jx( is+istride, js+jstride, ks+kstride);
                jy(i0+is, j0+js, k0+ks) += yee.jy( is+istride, js+jstride, ks+kstride);
                jz(i0+is, j0+js, k0+ks) += yee.jz( is+istride, js+jstride, ks+kstride);
                                                                                      
                rh(i0+is, j0+js, k0+ks) += yee.rho(is+istride, js+jstride, ks+kstride);
              }
            }
          }
        }
    //  }
    //}
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
    "/quick-fields-" +
    std::to_string(grid.comm.rank()) + 
    "_" +
    std::to_string(lap) +
    extension;
  std::cout << "QW: " << full_filename << std::endl;


  // open file and write
  File file(full_filename, H5F_ACC_TRUNC);
  file["Nx"] = arrs[0].Nx;
  file["Ny"] = arrs[0].Ny;
  file["Nz"] = arrs[0].Nz;

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

  return true;
}






//--------------------------------------------------
// explicit template class instantiations
//template class h5io::QuickWriter<1>;
template class h5io::QuickWriter<2>;
