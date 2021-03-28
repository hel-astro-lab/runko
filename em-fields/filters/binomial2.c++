#include "binomial2.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


/// single 2D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D 3-point binomial coefficients
  //real_short winv=1./16.,
  //        wtm=4.*winv, //middle
  //        wts=2.*winv, //side
  //        wtc=1.*winv; //corner
    
  // 2D 3-point binomial coefficients
  const real_short C2[3][3] = 
        { {1./16., 2./16., 1./16.},
          {2./16., 4./16., 2./16.},
          {1./16., 2./16., 1./16.} };

  auto& mesh = tile.get_yee();

  // using tmp as scratch arrays
  //
  // values are processed in the following order
  //
  // | 7 8 9 |
  // | 4 5 6 |
  // | 1 2 3 |
  // 
  // because this access pattern translates to continuous sweeps of
  // ...,1,2,3,....,4,5,6,....,7,8,9,.....
  // for the target array memory.

  //--------------------------------------------------
  // Jx

  const int halo = 2;

  const int imin = 0 - halo;
  const int imax = tile.mesh_lengths[0] + halo;

  const int jmin = 0 - halo;
  const int jmax = tile.mesh_lengths[1] + halo;

  // NOTE: using tmp as scratch arrays
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  for(int is=-1; is<=1; is++) 
  for(int js=-1; js<=1; js++) {
    real_short C = C2[is+1][js+1];

    for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,0) += mesh.jx(i+is, j+js, 0)*C;
        }
      }
  }
  swap(mesh.jx, tmp);
  
  //TODO: implement
  //assert(false);

}


/// single 3D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<3>::solve(
    fields::Tile<3>& tile)
{

  // 3D 3-point binomial coefficients
  const real_short C3[3][3][3] = 
        { { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} },
          { {2./64., 4./64., 2./64.}, {4./64., 8./64., 4./64.}, {2./64., 4./64., 2./64.} },
          { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} } };
  real_short C;

  //real_short winv  = 1./64., // normalization
  //           wtd  = 1.*winv, // diagnoal
  //           wtos = 2.*winv, // outer side
  //           wtis = 4.*winv, // inner side
  //           wt   = 8.*winv; // center
  //

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& mesh = tile.get_yee();
  const int halo = 2; 

  const int imin = 0 - halo;
  const int imax = tile.mesh_lengths[0] + halo;

  const int jmin = 0 - halo;
  const int jmax = tile.mesh_lengths[1] + halo;

  const int kmin = 0 - halo;
  const int kmax = tile.mesh_lengths[2] + halo;

  //--------------------------------------------------
  // Jx

  // make 3d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, int k, toolbox::Mesh<real_short, 3> &jj, toolbox::Mesh<real_short, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
    for(int js=-1; js<=1; js++) {
    for(int ks=-1; ks<=1; ks++) {
      tmp(i,j,k) += jj(i+is, j+js, k+ks) * C3[is+1][js+1][ks+1];
    }}}

  };

  // TODO: check that new 3x3x3 loop is equal to previous version
  tmp.clear();
  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jy, tmp);

  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jz, tmp);

  UniIter::sync();
  std::swap(mesh.jz, tmp);
  

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}


//template class fields::Binomial2<1>; // 1D
template class fields::Binomial2<2>; // 2D
template class fields::Binomial2<3>; // 3D
