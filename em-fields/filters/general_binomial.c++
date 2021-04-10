#include "general_binomial.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


template<>
void fields::General3p<2>::solve(
    fields::Tile<2>& tile)
{

  // 2D general coefficients
  const double winv=1./4.;                         //normalization
  const double wtm=winv*4.0*alpha*alpha,         //middle
               wts=winv*2.0*alpha*(1.0-alpha),   //side
               wtc=winv*(1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;


  const int k = 0;

  //--------------------------------------------------

  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    tmp(i,j,k) = 
        jj(i-1, j-1, k)*wtc + 
        jj(i  , j-1, k)*wts + 
        jj(i+1, j-1, k)*wtc + 
        jj(i-1, j  , k)*wts + 
        jj(i  , j  , k)*wtm + 
        jj(i+1, j  , k)*wts + 
        jj(i-1, j+1, k)*wtc + 
        jj(i  , j+1, k)*wts + 
        jj(i+1, j+1, k)*wtc;
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


//template class fields::General3p<1>; // 1D
template class fields::General3p<2>; // 2D
//template class fields::General3p<3>; // 3D
