#include <mpi4cpp/mpi.h>

#include "io/snapshots/fields.h"
#include "external/ezh5/src/ezh5.hpp"
#include "core/emf/tile.h"


using namespace mpi4cpp;
using ezh5::File;

template<>
inline void
  h5io::FieldsWriter<3>::read_tiles(corgi::Grid<3>& grid)
{

  // clear target arrays
  for(auto& arr: arrs) arr.clear();

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
  for(auto cid: grid.get_local_tiles()) {

    emf::Tile<3>& tile = [&]() -> emf::Tile<3>& {
      try {
        return dynamic_cast<emf::Tile<3>&>(grid.get_tile(cid));
      } catch(const std::bad_cast& ex) {
        throw std::runtime_error { std::format(
          "h5io::FieldsWriter::read_tiles assumes that all tiles are"
          "emf::Tiles or its descendants. Orginal exception: {}",
          ex.what()) };
      }
    }();

    const auto [Emds, Bmds, Jmds] = tile.view_EBJ_on_host();
    const auto [Nxuz, Nyuz, Nzuz] = tile.extents_wout_halo();

    const auto Nx = static_cast<int>(Nxuz);
    const auto Ny = static_cast<int>(Nyuz);
    const auto Nz = static_cast<int>(Nzuz);

    // starting location
    const int i0 = (Nx / stride) * std::get<0>(tile.index);
    const int j0 = (Ny / stride) * std::get<1>(tile.index);
    const int k0 = (Nz / stride) * std::get<2>(tile.index);

    const int nxt = Nx / stride;
    const int nyt = Ny / stride;
    const int nzt = Nz / stride;

    // copy tile patch by stride hopping; either downsample or average

    // field quantities; just downsample by hopping with stride
    for(int ks = 0; ks < nzt; ks++)
      for(int js = 0; js < nyt; js++)
        for(int is = 0; is < nxt; is++) {
          ex(i0 + is, j0 + js, k0 + ks) =
            Emds[is * stride, js * stride, ks * stride][0];
          ey(i0 + is, j0 + js, k0 + ks) =
            Emds[is * stride, js * stride, ks * stride][1];
          ez(i0 + is, j0 + js, k0 + ks) =
            Emds[is * stride, js * stride, ks * stride][2];

          bx(i0 + is, j0 + js, k0 + ks) =
            Bmds[is * stride, js * stride, ks * stride][0];
          by(i0 + is, j0 + js, k0 + ks) =
            Bmds[is * stride, js * stride, ks * stride][1];
          bz(i0 + is, j0 + js, k0 + ks) =
            Bmds[is * stride, js * stride, ks * stride][2];
        }

    // densities; these quantities we average over the volume
    for(int ks = 0; ks < nzt; ks++)
      for(int kstride = 0; kstride < stride; kstride++)
        for(int js = 0; js < nyt; js++)
          for(int jstride = 0; jstride < stride; jstride++)
            for(int is = 0; is < nxt; is++)
              for(int istride = 0; istride < stride; istride++) {
                jx(i0 + is, j0 + js, k0 + ks) += Jmds
                  [is * stride + istride, js * stride + jstride, ks * stride + kstride]
                  [0];
                jy(i0 + is, j0 + js, k0 + ks) += Jmds
                  [is * stride + istride, js * stride + jstride, ks * stride + kstride]
                  [1];
                jz(i0 + is, j0 + js, k0 + ks) += Jmds
                  [is * stride + istride, js * stride + jstride, ks * stride + kstride]
                  [2];
                rh(i0 + is, j0 + js, k0 + ks) += 0; /* gs.rho(
                  is * stride + istride,
                  js * stride + jstride,
                  ks * stride + kstride); */
              }

  }  // tiles
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
//template class h5io::FieldsWriter<2>;
template class h5io::FieldsWriter<3>;
