#include "general_binomial.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


void sweep_in_x(
        toolbox::Mesh<float_m,3>& arr,
        int lenx,
        int leny,
        int lenz,
        int sub_cycle
        ) {

  int imin, imax, jmin, jmax, kmin, kmax;
  float_m tmp1, tmp2;
    
  for(int npass=1; npass<=sub_cycle; npass++){
    imin = -3 + npass;
    imax = lenx + 3 - 1 - npass;

    jmin = -3 + npass;
    jmax = leny + 3 - 1 - npass;

    kmin = -3 + npass;
    kmax = lenz + 3 - 1 - npass;

    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {

        tmp2 = arr(imin-1, j, k);
        tmp2 = arr(imin-1, j, k);
        for(int i=imin; i<imax-1; i=i+2) {
        //for i = imin, imax - 1, 2 //mod 2 imax-imin
        //for i = imin, imax, 2     //non mod2 

          tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k);

          arr(i - 1, j, k) = tmp2;
          tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i + 1, j, k) + 0.25 * arr(i + 2, j, k);
          arr(i, j, k) = tmp1;
          }

        // extra manual step for uneven grids
        tmp1 = 0.25 * arr(imax - 1, j, k) + 0.5 * arr(imax, j, k) + 0.25 * arr(imax + 1, j, k);
        arr(imax - 1, j, k) = tmp2;

        arr(imax, j, k) = tmp1;
      }//j
    }//k
  } // end of npass
}



/// Single 3D 2nd order 3-point Binomial filter
// NOTE: optimized to slide in x/y/z dimensions to enable vectorization
// TODO: not fully implemented but frame is presented here
template<>
void emf::OptBinomial2<3>::solve(
    emf::Tile<3>& tile)
{
  auto& mesh = tile.get_grids();

  sweep_in_x(mesh.jx, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);
  sweep_in_x(mesh.jy, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);
  sweep_in_x(mesh.jz, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);

  //sweep_in_y(mesh.jx, 1);
  //sweep_in_y(mesh.jy, 1);
  //sweep_in_y(mesh.jz, 1);

  //sweep_in_z(mesh.jx, 1);
  //sweep_in_z(mesh.jy, 1);
  //sweep_in_z(mesh.jz, 1);

}


//template class emf::OptBinomial2<1>; // 1D
  template class emf::OptBinomial2<2>; // 2D
//template class emf::OptBinomial2<3>; // 3D
