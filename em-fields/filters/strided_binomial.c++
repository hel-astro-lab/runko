#include "strided_binomial.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


template<>
void fields::General3pStrided<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  const double winv=1./4.;                         //normalization
  const double wtm=winv * 4.0*alpha*alpha,         //middle
               wts=winv * 2.0*alpha*(1.0-alpha),   //side
               wtc=winv * (1.0-alpha)*(1.0-alpha); //corner
    
  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;

  const int k = 0;
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    for(int istr = 1; istr < stride; istr++){
      tmp(i,j,k) = 
          jj(i-istr, j-istr, k)*wtc + 
          jj(i     , j-istr, k)*wts + 
          jj(i+istr, j-istr, k)*wtc + 

          jj(i-istr, j     , k)*wts + 
          jj(i     , j     , k)*wtm + 
          jj(i+istr, j     , k)*wts + 

          jj(i-istr, j+istr, k)*wtc + 
          jj(i     , j+istr, k)*wts + 
          jj(i+istr, j+istr, k)*wtc;
      }
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jy, tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jz, tmp);
 
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
void fields::Binomial2Strided2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  const double wn=1./16.0/16.0;  //normalization
    
  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;

  const int k = 0;
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    tmp(i,j,k) = 
        jj(i-3, j-3, 0)*1.0*wn +
        jj(i-2, j-3, 0)*2.0*wn +
        jj(i-1, j-3, 0)*3.0*wn +
        jj(i+0, j-3, 0)*4.0*wn +
        jj(i+1, j-3, 0)*3.0*wn +
        jj(i+2, j-3, 0)*2.0*wn +
        jj(i+3, j-3, 0)*1.0*wn +
        jj(i-3, j-2, 0)*2.0*wn +
        jj(i-2, j-2, 0)*4.0*wn +
        jj(i-1, j-2, 0)*6.0*wn +
        jj(i+0, j-2, 0)*8.0*wn +
        jj(i+1, j-2, 0)*6.0*wn +
        jj(i+2, j-2, 0)*4.0*wn +
        jj(i+3, j-2, 0)*2.0*wn +
        jj(i-3, j-1, 0)*3.0*wn +
        jj(i-2, j-1, 0)*6.0*wn +
        jj(i-1, j-1, 0)*9.0*wn +
        jj(i+0, j-1, 0)*12.*wn +
        jj(i+1, j-1, 0)*9.0*wn +
        jj(i+2, j-1, 0)*6.0*wn +
        jj(i+3, j-1, 0)*3.0*wn +
        jj(i-3, j+0, 0)*4.0*wn +
        jj(i-2, j+0, 0)*8.0*wn +
        jj(i-1, j+0, 0)*12.*wn +
        jj(i+0, j+0, 0)*16.*wn +
        jj(i+1, j+0, 0)*12.*wn +
        jj(i+2, j+0, 0)*8.0*wn +
        jj(i+3, j+0, 0)*4.0*wn +
        jj(i-3, j+1, 0)*3.0*wn +
        jj(i-2, j+1, 0)*6.0*wn +
        jj(i-1, j+1, 0)*9.0*wn +
        jj(i+0, j+1, 0)*12.*wn +
        jj(i+1, j+1, 0)*9.0*wn +
        jj(i+2, j+1, 0)*6.0*wn +
        jj(i+3, j+1, 0)*3.0*wn +
        jj(i-3, j+2, 0)*2.0*wn +
        jj(i-2, j+2, 0)*4.0*wn +
        jj(i-1, j+2, 0)*6.0*wn +
        jj(i+0, j+2, 0)*8.0*wn +
        jj(i+1, j+2, 0)*6.0*wn +
        jj(i+2, j+2, 0)*4.0*wn +
        jj(i+3, j+2, 0)*2.0*wn +
        jj(i-3, j+3, 0)*1.0*wn +
        jj(i-2, j+3, 0)*2.0*wn +
        jj(i-1, j+3, 0)*3.0*wn +
        jj(i+0, j+3, 0)*4.0*wn +
        jj(i+1, j+3, 0)*3.0*wn +
        jj(i+2, j+3, 0)*2.0*wn +
        jj(i+3, j+3, 0)*1.0*wn;
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jy, tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate2D(fun, 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jz, tmp);
 
  UniIter::sync();
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}



//template class fields::General3pStrided<1>; // 1D
template class fields::General3pStrided<2>; // 2D
//template class fields::General3pStrided<3>; // 3D
  
//template class fields::Binomial2Strided2<1>; // 1D
template class fields::Binomial2Strided2<2>; // 2D
//template class fields::Binomial2Strided2<3>; // 3D
