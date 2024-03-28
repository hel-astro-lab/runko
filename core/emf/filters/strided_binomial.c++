
#include <cmath>

#include "core/emf/filters/strided_binomial.h"
#include "external/iter/devcall.h"
#include "external/iter/iter.h"
#include "external/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


template<>
void emf::General3pStrided<2>::solve(
    emf::Tile<2>& tile)
{
  // 2D general coefficients
  const double winv=1./4.;                         //normalization
  const double wtm=winv * 4.0*alpha*alpha,         //middle
               wts=winv * 2.0*alpha*(1.0-alpha),   //side
               wtc=winv * (1.0-alpha)*(1.0-alpha); //corner
    
  auto& mesh = tile.get_grids();

  const int H = 0; 
  const int k = 0;
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float, 3> &jj, 
                   toolbox::Mesh<float, 3> &tmp)
  {
    for(int istr = 1; istr < stride; istr++){
      tmp(i-H,j-H,k) = 
          jj(i-istr-H, j-istr-H, k)*wtc + 
          jj(i     -H, j-istr-H, k)*wts + 
          jj(i+istr-H, j-istr-H, k)*wtc + 

          jj(i-istr-H, j     -H, k)*wts + 
          jj(i     -H, j     -H, k)*wtm + 
          jj(i+istr-H, j     -H, k)*wts + 

          jj(i-istr-H, j+istr-H, k)*wtc + 
          jj(i     -H, j+istr-H, k)*wts + 
          jj(i+istr-H, j+istr-H, k)*wtc;
      }
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jy, tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jz, tmp);
 
  UniIter::sync();
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}



/// outwinded 2-strided 3-point binomial
//  Optimal for 3-point halo regions; 
//  however on some platforms looping over strides (i.e. general filter)
//  might be faster because of prefetcher behavior
//  
//  NOTE: Simple timing shows that this is 2x slower than general form.
//
// coefficients are:
// [ 1.  2.  3.  4.  3.  2.  1. ]
// [ 2.  4.  6.  8.  6.  4.  2. ]
// [ 3.  6.  9. 12.  9.  6.  3. ]
// [ 4.  8. 12. 16. 12.  8.  4. ]
// [ 3.  6.  9. 12.  9.  6.  3. ]
// [ 2.  4.  6.  8.  6.  4.  2. ]
// [ 1.  2.  3.  4.  3.  2.  1. ]
template<>
void emf::Binomial2Strided2<2>::solve(
    emf::Tile<2>& tile)
{
  // 2D general coefficients
  const double wn=1./16.0/16.0;  //normalization
    
  auto& mesh = tile.get_grids();
  const int H = 0; 
  const int k = 0;
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float, 3> &jj, 
                   toolbox::Mesh<float, 3> &tmp)
  {
    tmp(i-H,j-H,k) = 
        jj(i-3-H, j-3-H, 0)*1.0*wn +
        jj(i-2-H, j-3-H, 0)*2.0*wn +
        jj(i-1-H, j-3-H, 0)*3.0*wn +
        jj(i+0-H, j-3-H, 0)*4.0*wn +
        jj(i+1-H, j-3-H, 0)*3.0*wn +
        jj(i+2-H, j-3-H, 0)*2.0*wn +
        jj(i+3-H, j-3-H, 0)*1.0*wn +
        jj(i-3-H, j-2-H, 0)*2.0*wn +
        jj(i-2-H, j-2-H, 0)*4.0*wn +
        jj(i-1-H, j-2-H, 0)*6.0*wn +
        jj(i+0-H, j-2-H, 0)*8.0*wn +
        jj(i+1-H, j-2-H, 0)*6.0*wn +
        jj(i+2-H, j-2-H, 0)*4.0*wn +
        jj(i+3-H, j-2-H, 0)*2.0*wn +
        jj(i-3-H, j-1-H, 0)*3.0*wn +
        jj(i-2-H, j-1-H, 0)*6.0*wn +
        jj(i-1-H, j-1-H, 0)*9.0*wn +
        jj(i+0-H, j-1-H, 0)*12.*wn +
        jj(i+1-H, j-1-H, 0)*9.0*wn +
        jj(i+2-H, j-1-H, 0)*6.0*wn +
        jj(i+3-H, j-1-H, 0)*3.0*wn +
        jj(i-3-H, j+0-H, 0)*4.0*wn +
        jj(i-2-H, j+0-H, 0)*8.0*wn +
        jj(i-1-H, j+0-H, 0)*12.*wn +
        jj(i+0-H, j+0-H, 0)*16.*wn +
        jj(i+1-H, j+0-H, 0)*12.*wn +
        jj(i+2-H, j+0-H, 0)*8.0*wn +
        jj(i+3-H, j+0-H, 0)*4.0*wn +
        jj(i-3-H, j+1-H, 0)*3.0*wn +
        jj(i-2-H, j+1-H, 0)*6.0*wn +
        jj(i-1-H, j+1-H, 0)*9.0*wn +
        jj(i+0-H, j+1-H, 0)*12.*wn +
        jj(i+1-H, j+1-H, 0)*9.0*wn +
        jj(i+2-H, j+1-H, 0)*6.0*wn +
        jj(i+3-H, j+1-H, 0)*3.0*wn +
        jj(i-3-H, j+2-H, 0)*2.0*wn +
        jj(i-2-H, j+2-H, 0)*4.0*wn +
        jj(i-1-H, j+2-H, 0)*6.0*wn +
        jj(i+0-H, j+2-H, 0)*8.0*wn +
        jj(i+1-H, j+2-H, 0)*6.0*wn +
        jj(i+2-H, j+2-H, 0)*4.0*wn +
        jj(i+3-H, j+2-H, 0)*2.0*wn +
        jj(i-3-H, j+3-H, 0)*1.0*wn +
        jj(i-2-H, j+3-H, 0)*2.0*wn +
        jj(i-1-H, j+3-H, 0)*3.0*wn +
        jj(i+0-H, j+3-H, 0)*4.0*wn +
        jj(i+1-H, j+3-H, 0)*3.0*wn +
        jj(i+2-H, j+3-H, 0)*2.0*wn +
        jj(i+3-H, j+3-H, 0)*1.0*wn;
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jy, tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        mesh.jz, tmp);
 
  UniIter::sync();
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}



//template class emf::General3pStrided<1>; // 1D
template class emf::General3pStrided<2>; // 2D
//template class emf::General3pStrided<3>; // 3D
  
//template class emf::Binomial2Strided2<1>; // 1D
template class emf::Binomial2Strided2<2>; // 2D
//template class emf::Binomial2Strided2<3>; // 3D
