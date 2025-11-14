
#include <cmath>

#include "core/emf/filters/binomial2.h"
#include "external/iter/devcall.h"
#include "external/iter/iter.h"
#include "external/iter/allocator.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


/// single 1D 2nd order 3-point binomial filter 
template<>
void emf::Binomial2<1>::solve(
    emf::Tile<1>& tile)
{
    
  // 1D 3-point binomial coefficients
  const float C1[3] = {1./4., 2./4., 1./4.};

  auto& mesh = tile.get_grids();

  // halo width
  const int H = 2; 

  // NOTE: using tmp as scratch arrays

  // make 1d loop with shared memory 
  // NOTE: shifted with -H to iterate over halos
  // NOTE: similarly, limits are expanded by 2*H
  auto fun = 
  [=] DEVCALLABLE (int i,  
                   toolbox::Mesh<float, 3> &jj, 
                   toolbox::Mesh<float, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
      tmp(i-H,0,0) += jj(i+is-H, 0, 0)*C1[is+1];
      //tmp(i) += jj(i+is)*C1[is+1]; // NOTE: raw 1d indexing does not work
    }
  };
    
  //--------------------------------------------------
  // Jx
  tmp.clear();
  UniIter::iterate(fun, 
        tile.mesh_lengths[0] + 2*H, 
        mesh.jx, 
        tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  UniIter::iterate(fun, 
        tile.mesh_lengths[0] + 2*H, 
        mesh.jy, 
        tmp);
 
  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  UniIter::iterate(fun, 
        tile.mesh_lengths[0] + 2*H, 
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
void emf::Binomial2<2>::solve(
    emf::Tile<2>& tile)
{
    
  // 2D 3-point binomial coefficients
  const float C2[3][3] = 
        { {1./16., 2./16., 1./16.},
          {2./16., 4./16., 2./16.},
          {1./16., 2./16., 1./16.} };

  auto& mesh = tile.get_grids();

  const int H = 2; 

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


  // NOTE: using tmp as scratch arrays
    
  // make 2d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, 
                   toolbox::Mesh<float, 3> &jj, 
                   toolbox::Mesh<float, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
    for(int js=-1; js<=1; js++) {
      tmp(i-H,j-H,0) += jj(i+is-H, j+js-H, 0)*C2[is+1][js+1];
    }}
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


/// single 3D 2nd order 3-point binomial filter 
template<>
void emf::Binomial2<3>::solve(
    emf::Tile<3>& tile)
{
#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif


  // 3D 3-point binomial coefficients
  const float C3[3][3][3] = 
        { { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} },
          { {2./64., 4./64., 2./64.}, {4./64., 8./64., 4./64.}, {2./64., 4./64., 2./64.} },
          { {1./64., 2./64., 1./64.}, {2./64., 4./64., 2./64.}, {1./64., 2./64., 1./64.} } };

  auto& mesh = tile.get_grids();
  const int H = 2; 


  // make 3d loop with shared memory 
  auto fun = 
  [=] DEVCALLABLE (int i, int j, int k, toolbox::Mesh<float, 3> &jj, toolbox::Mesh<float, 3> &tmp)
  {
    for(int is=-1; is<=1; is++) {
    for(int js=-1; js<=1; js++) {
    for(int ks=-1; ks<=1; ks++) {
      tmp(i-H,j-H,k-H) += jj(i+is-H, j+js-H, k+ks-H) * C3[is+1][js+1][ks+1];
    }}}

  };

  tmp.clear();
  UniIter::iterate3D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        tile.mesh_lengths[2] + 2*H,
        mesh.jx, tmp);
 
  UniIter::sync();
  std::swap(mesh.jx, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        tile.mesh_lengths[2] + 2*H,
        mesh.jy, tmp);

  UniIter::sync();
  std::swap(mesh.jy, tmp);

  //--------------------------------------------------
  tmp.clear();
  UniIter::iterate3D(fun, 
        tile.mesh_lengths[0] + 2*H, 
        tile.mesh_lengths[1] + 2*H,
        tile.mesh_lengths[2] + 2*H,
        mesh.jz, tmp);

  UniIter::sync();
  std::swap(mesh.jz, tmp);
  

  //--------------------------------------------------
#ifdef GPU
  nvtxRangePop();
#endif
}


template class emf::Binomial2<1>; // 1D
template class emf::Binomial2<2>; // 2D
template class emf::Binomial2<3>; // 3D
