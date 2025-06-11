
#include <cmath>

#include "core/emf/filters/general_binomial.h"
#include "external/iter/iter.h"

template<>
void emf::General3p<2>::solve(
    emf::Tile<2>& tile)
{

  // 2D general coefficients
  const double winv=1./4.;                         //normalization
  const double wtm=winv*4.0*alpha*alpha,         //middle
               wts=winv*2.0*alpha*(1.0-alpha),   //side
               wtc=winv*(1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_grids();
  const int H = 2; 
  const int k = 0;

  //--------------------------------------------------

  // make 2d loop with shared memory 
  auto fun = 
  [=]  (int i, int j, 
                   toolbox::Mesh<float, 3> &jj, 
                   toolbox::Mesh<float, 3> &tmp)
  {
    tmp(i-H,j-H,k) = 
        jj(i-1-H, j-1-H, k)*wtc + 
        jj(i  -H, j-1-H, k)*wts + 
        jj(i+1-H, j-1-H, k)*wtc + 
        jj(i-1-H, j  -H, k)*wts + 
        jj(i  -H, j  -H, k)*wtm + 
        jj(i+1-H, j  -H, k)*wts + 
        jj(i-1-H, j+1-H, k)*wtc + 
        jj(i  -H, j+1-H, k)*wts + 
        jj(i+1-H, j+1-H, k)*wtc;
  };

  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jx, tmp);
 
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jy, tmp);
 
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jz, tmp);
 
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
}


//template class emf::General3p<1>; // 1D
template class emf::General3p<2>; // 2D
//template class emf::General3p<3>; // 3D
