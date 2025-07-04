
#include <cmath>

#include "core/emf/filters/compensator.h"
#include "external/iter/iter.h"

/// 2D 3-point compensator filter from Birdsall & Langdon
//
// Filter is (20,-1,-1) with normalization 1/12
//
//  NOTE: Main difference to a binomial compensator are
//  the suppressed corners
//
template<>
void emf::Compensator2<2>::solve(
    emf::Tile<2>& tile)
{
  // 2D general coefficients
  const double winv=1./12.; //normalization
  const double wtm=20.0*winv, //middle M
               wts=-1.0*winv, //side   S
               wtc=-1.0*winv; //corner K


  auto& mesh = tile.get_grids();
  const int H = 2; 
  const int k = 0;

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



//template class emf::Compensator2<1>; // 1D
  template class emf::Compensator2<2>; // 2D
//template class emf::Compensator2<3>; // 3D
