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
void fields::Binomial2<1>::solve(
    fields::Tile<1>& tile)
{
    
  // 2D 3-point binomial coefficients
  const float_m C1[3] = {1./4., 2./4., 1./4.};

  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int imax = tile.mesh_lengths[0] + halo;


  //--------------------------------------------------
  // Jx

  // NOTE: using tmp as scratch arrays
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i,  
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
      tmp(i,0,0) += jj(i+is, 0, 0)*C1[is+1];
    }
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        mesh.jx, 
        tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        mesh.jy, 
        tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        mesh.jz, 
        tmp);
 
  UniIter::sync();
  std::swap(mesh.jz, tmp);

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}


/// single 2D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<2>::solve(
    fields::Tile<2>& tile)
{
    
  // 2D 3-point binomial coefficients
  const float_m C2[3][3] = 
        { {1./16., 2./16., 1./16.},
          {2./16., 4./16., 2./16.},
          {1./16., 2./16., 1./16.} };

  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;


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

  // NOTE: using tmp as scratch arrays
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float_m, 3> &jj, 
                   toolbox::Mesh<float_m, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
    for(int js=-1; js<=1; js++) {
      tmp(i,j,0) += jj(i+is, j+js, 0)*C2[is+1][js+1];
    }}
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


/// single 3D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<3>::solve(
    fields::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif


  // 3D 3-point binomial coefficients
  const float_m C3[3][3][3] = 
        { { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} },
          { {2./64., 4./64., 2./64.}, {4./64., 8./64., 4./64.}, {2./64., 4./64., 2./64.} },
          { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} } };

  auto& mesh = tile.get_yee();

  const int halo = 2; 

  const int imin = 0 - halo;
  const int jmin = 0 - halo;
  const int kmin = 0 - halo;

  const int imax = tile.mesh_lengths[0] + halo;
  const int jmax = tile.mesh_lengths[1] + halo;
  const int kmax = tile.mesh_lengths[2] + halo;

  //--------------------------------------------------
  // Jx

  // make 3d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, int k, toolbox::Mesh<float_m, 3> &jj, toolbox::Mesh<float_m, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
    for(int js=-1; js<=1; js++) {
    for(int ks=-1; ks<=1; ks++) {
      tmp(i,j,k) += jj(i+is, j+js, k+ks) * C3[is+1][js+1][ks+1];
    }}}

  };

  // TODO: check that new 3x3x3 loop is equal to previous version
  tmp.clear();
  UniIter::iterate3D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[2]),
        mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[2]),
        mesh.jy, tmp);

  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, 
        static_cast<int>(tile.mesh_lengths[0]), 
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[2]),
        mesh.jz, tmp);

  UniIter::sync();
  std::swap(mesh.jz, tmp);
  

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}


template class fields::Binomial2<1>; // 1D
template class fields::Binomial2<2>; // 2D
template class fields::Binomial2<3>; // 3D
