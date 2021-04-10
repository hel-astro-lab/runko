#include "compensator.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


/// 2D 3-point compensator filter from Birdsall & Langdon
//
// Filter is (20,-1,-1) with normalization 1/12
//
//  NOTE: Main difference to a binomial compensator are
//  the suppressed corners
//
template<>
void fields::Compensator2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  const double winv=1./12.; //normalization
  const double wtm=20.0*winv, //middle M
               wts=-1.0*winv, //side   S
               wtc=-1.0*winv; //corner K


  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;


  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    tmp(i,j,0) = 
      mesh.jx(i-1, j-1, 0)*wtc + 
      mesh.jx(i  , j-1, 0)*wts + 
      mesh.jx(i+1, j-1, 0)*wtc + 

      mesh.jx(i-1, j  , 0)*wts + 
      mesh.jx(i  , j  , 0)*wtm + 
      mesh.jx(i+1, j  , 0)*wts + 

      mesh.jx(i-1, j+1, 0)*wtc + 
      mesh.jx(i  , j+1, 0)*wts + 
      mesh.jx(i+1, j+1, 0)*wtc;
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        mesh.jy, tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        mesh.jz, tmp);
 
  UniIter::sync();
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}



//template class fields::Compensator2<1>; // 1D
  template class fields::Compensator2<2>; // 2D
//template class fields::Compensator2<3>; // 3D
