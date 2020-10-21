#include "digital.h"

#include <cmath>
#include "../../tools/iter/devcall.h"
#include "../../tools/iter/iter.h"
#include "../../tools/iter/allocator.h"

#include <nvtx3/nvToolsExt.h> 

/// single 2D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D 3-point binomial coefficients
  real_short winv=1./16.,
          wtm=4.*winv, //middle
          wts=2.*winv, //side
          wtc=1.*winv; //corner

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
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jx(i-1, j-1, k)*wtc + 
      mesh.jx(i  , j-1, k)*wts + 
      mesh.jx(i+1, j-1, k)*wtc + 

      mesh.jx(i-1, j  , k)*wts + 
      mesh.jx(i  , j  , k)*wtm + 
      mesh.jx(i+1, j  , k)*wts + 

      mesh.jx(i-1, j+1, k)*wtc + 
      mesh.jx(i  , j+1, k)*wts + 
      mesh.jx(i+1, j+1, k)*wtc;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jy(i-1, j-1, k)*wtc + 
      mesh.jy(i  , j-1, k)*wts + 
      mesh.jy(i+1, j-1, k)*wtc + 

      mesh.jy(i-1, j  , k)*wts + 
      mesh.jy(i  , j  , k)*wtm + 
      mesh.jy(i+1, j  , k)*wts + 

      mesh.jy(i-1, j+1, k)*wtc + 
      mesh.jy(i  , j+1, k)*wts + 
      mesh.jy(i+1, j+1, k)*wtc;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jz(i-1, j-1, k)*wtc + 
      mesh.jz(i  , j-1, k)*wts + 
      mesh.jz(i+1, j-1, k)*wtc + 

      mesh.jz(i-1, j  , k)*wts + 
      mesh.jz(i  , j  , k)*wtm + 
      mesh.jz(i+1, j  , k)*wts + 

      mesh.jz(i-1, j+1, k)*wtc + 
      mesh.jz(i  , j+1, k)*wts + 
      mesh.jz(i+1, j+1, k)*wtc;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}

//#define GPU
#ifdef GPU

    template<class F, class... Args>
    __global__ void iterate3dShared(F fun, int xMax, int yMax, int zMax, toolbox::Mesh<real_short, 3> *jj, toolbox::Mesh<real_short, 3> *tmp)
    {
        real_short winv  = 1./64., // normalization
             wtd  = 1.*winv, // diagnoal
             wtos = 2.*winv, // outer side
             wtis = 4.*winv, // inner side
             wt   = 8.*winv; // center

        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int idy = blockIdx.y * blockDim.y + threadIdx.y;
        int idz = blockIdx.z * blockDim.z + threadIdx.z;

        int locIdx = threadIdx.x;
        int locIdy = threadIdx.y;
        int locIdz = threadIdx.z;

        __shared__ real_short local[10][10][10];


          for (int x = threadIdx.x; x < (blockDim.x+2); x += blockDim.x)
          for (int y = threadIdx.y; y < (blockDim.y+2); y += blockDim.y)
          for (int z = threadIdx.z; z < (blockDim.z+2); z += blockDim.z)
        {
          int targetX = (blockIdx.x * blockDim.x) + (x-1);
          int targetY = (blockIdx.y * blockDim.y) + (y-1);
          int targetZ = (blockIdx.z * blockDim.z) + (z-1);
          local[x][y][z] = (*jj)(
            targetX, 
            targetY, 
            targetZ);
        }
        
      
        __syncthreads();
        
        if(idx >= xMax) return;
        if(idy >= yMax) return;
        if(idz >= zMax) return;

        locIdx++;
        locIdy++;
        locIdz++;

        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = blockIdx.y * blockDim.y + threadIdx.y;
        int k = blockIdx.z * blockDim.z + threadIdx.z;

    		(*tmp)(idx,idy,idz) =
          local[locIdx-1][locIdy-1][locIdz-1]*wtd  + 
          local[locIdx  ][locIdy-1][locIdz-1]*wtos + 
          local[locIdx+1][locIdy-1][locIdz-1]*wtd  +
             
          local[locIdx-1][locIdy  ][locIdz-1]*wtos +
          local[locIdx  ][locIdy  ][locIdz-1]*wtis +
          local[locIdx+1][locIdy  ][locIdz-1]*wtos +
             
          local[locIdx-1][locIdy+1][locIdz-1]*wtd  + 
          local[locIdx  ][locIdy+1][locIdz-1]*wtos +
          local[locIdx+1][locIdy+1][locIdz-1]*wtd  +
             
          local[locIdx-1][locIdy-1][locIdz  ]*wtos +
          local[locIdx  ][locIdy-1][locIdz  ]*wtis +
          local[locIdx+1][locIdy-1][locIdz  ]*wtos +
             
          local[locIdx-1][locIdy  ][locIdz  ]*wtis +
          local[locIdx  ][locIdy  ][locIdz  ]*wt   +
          local[locIdx+1][locIdy  ][locIdz  ]*wtis +
             
          local[locIdx-1][locIdy+1][locIdz  ]*wtos +
          local[locIdx  ][locIdy+1][locIdz  ]*wtis +
          local[locIdx+1][locIdy+1][locIdz  ]*wtos +
             
          local[locIdx-1][locIdy-1][locIdz+1]*wtd  +
          local[locIdx  ][locIdy-1][locIdz+1]*wtos +
          local[locIdx+1][locIdy-1][locIdz+1]*wtd  +
             
          local[locIdx-1][locIdy  ][locIdz+1]*wtos +
          local[locIdx  ][locIdy  ][locIdz+1]*wtis +
          local[locIdx+1][locIdy  ][locIdz+1]*wtos +
             
          local[locIdx-1][locIdy+1][locIdz+1]*wtd  +
          local[locIdx  ][locIdy+1][locIdz+1]*wtos +
          local[locIdx+1][locIdy+1][locIdz+1]*wtd;
          
    }
#endif

/// single 3D 2nd order 3-point binomial filter 
template<>
void fields::Binomial2<3>::solve(
    fields::Tile<3>& tile)
{
  // 3D 3-point binomial coefficients
  real_short winv  = 1./64., // normalization
             wtd  = 1.*winv, // diagnoal
             wtos = 2.*winv, // outer side
             wtis = 4.*winv, // inner side
             wt   = 8.*winv; // center
nvtxRangePush(__PRETTY_FUNCTION__);

  auto& mesh = tile.get_yee();

  // using tmp as scratch arrays
  //
  // TODO: can be optimized by changing to 3x loops and using 3x3 array; 
  //       this should vectorize easier.

  //--------------------------------------------------
  // Jx
  DEV_REGISTER

// make with shared memory 
    auto fun = 
  [=] DEVCALLABLE (int i, int j, int k, toolbox::Mesh<real_short, 3> &jj, toolbox::Mesh<real_short, 3> &tmp)
  {
    		tmp(i,j,k) =
      jj(i-1, j-1, k-1)*wtd  + 
      jj(i  , j-1, k-1)*wtos + 
      jj(i+1, j-1, k-1)*wtd  +
             
      jj(i-1, j  , k-1)*wtos +
      jj(i  , j  , k-1)*wtis +
      jj(i+1, j  , k-1)*wtos +
             
      jj(i-1, j+1, k-1)*wtd  + 
      jj(i  , j+1, k-1)*wtos +
      jj(i+1, j+1, k-1)*wtd  +
             
      jj(i-1, j-1, k  )*wtos +
      jj(i  , j-1, k  )*wtis +
      jj(i+1, j-1, k  )*wtos +
             
      jj(i-1, j  , k  )*wtis +
      jj(i  , j  , k  )*wt   +
      jj(i+1, j  , k  )*wtis +
             
      jj(i-1, j+1, k  )*wtos +
      jj(i  , j+1, k  )*wtis +
      jj(i+1, j+1, k  )*wtos +
             
      jj(i-1, j-1, k+1)*wtd  +
      jj(i  , j-1, k+1)*wtos +
      jj(i+1, j-1, k+1)*wtd  +
             
      jj(i-1, j  , k+1)*wtos +
      jj(i  , j  , k+1)*wtis +
      jj(i+1, j  , k+1)*wtos +
             
      jj(i-1, j+1, k+1)*wtd  +
      jj(i  , j+1, k+1)*wtos +
      jj(i+1, j+1, k+1)*wtd;
  };

/*
  auto gridArgs = UniIter::UniIterCU::getGrid3({
        static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0])},{8,8,4});

  getErrorCuda(((iterate3dShared<<<std::get<1>(gridArgs), std::get<0>(gridArgs)>>>(fun, 
        static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), 
        UniIter::UniIterCU::getDevPtr(mesh.jx), UniIter::UniIterCU::getDevPtr(tmp)))));
  UniIter::sync();
  //mesh.jx = tmp; // then copy from scratch to original arrays
  std::swap(mesh.jx, tmp);

  getErrorCuda(((iterate3dShared<<<std::get<1>(gridArgs), std::get<0>(gridArgs)>>>(fun, 
        static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), 
        UniIter::UniIterCU::getDevPtr(mesh.jy), UniIter::UniIterCU::getDevPtr(tmp)))));


  UniIter::sync();
  //mesh.jy = tmp; // then copy from scratch to original arrays
  std::swap(mesh.jy, tmp);


  getErrorCuda(((iterate3dShared<<<std::get<1>(gridArgs), std::get<0>(gridArgs)>>>(fun, 
        static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), 
        UniIter::UniIterCU::getDevPtr(mesh.jz), UniIter::UniIterCU::getDevPtr(tmp)))));


  UniIter::sync();
  //mesh.jz = tmp; // then copy from scratch to original arrays // calls the copy constructor and the data from tmp is actually copied back, a possibly better solution is to just swap the pointers
  std::swap(mesh.jz, tmp);

*/


  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jx, tmp);

  UniIter::sync();
  //mesh.jx = tmp; // then copy from scratch to original arrays
  std::swap(mesh.jx, tmp);

  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jy, tmp);

  UniIter::sync();
  //mesh.jy = tmp; // then copy from scratch to original arrays
  std::swap(mesh.jy, tmp);


  UniIter::iterate3D(fun, static_cast<int>(tile.mesh_lengths[2]),
        static_cast<int>(tile.mesh_lengths[1]),
        static_cast<int>(tile.mesh_lengths[0]), mesh.jz, tmp);

  UniIter::sync();
  //mesh.jz = tmp; // then copy from scratch to original arrays // calls the copy constructor and the data from tmp is actually copied back, a possibly better solution is to just swap the pointers
  std::swap(mesh.jz, tmp);
  

nvtxRangePop();

}



/// single 2D 3-point general filter pass
template<>
void fields::General3p<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  double winv=1./4.;                         //normalization
  double wtm=winv * 4.0*alpha*alpha,         //middle
         wts=winv * 2.0*alpha*(1.0-alpha),   //side
         wtc=winv * (1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_yee();


  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jx(i-1, j-1, k)*wtc + 
      mesh.jx(i  , j-1, k)*wts + 
      mesh.jx(i+1, j-1, k)*wtc + 

      mesh.jx(i-1, j  , k)*wts + 
      mesh.jx(i  , j  , k)*wtm + 
      mesh.jx(i+1, j  , k)*wts + 

      mesh.jx(i-1, j+1, k)*wtc + 
      mesh.jx(i  , j+1, k)*wts + 
      mesh.jx(i+1, j+1, k)*wtc;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jy(i-1, j-1, k)*wtc + 
      mesh.jy(i  , j-1, k)*wts + 
      mesh.jy(i+1, j-1, k)*wtc + 

      mesh.jy(i-1, j  , k)*wts + 
      mesh.jy(i  , j  , k)*wtm + 
      mesh.jy(i+1, j  , k)*wts + 

      mesh.jy(i-1, j+1, k)*wtc + 
      mesh.jy(i  , j+1, k)*wts + 
      mesh.jy(i+1, j+1, k)*wtc;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jz(i-1, j-1, k)*wtc + 
      mesh.jz(i  , j-1, k)*wts + 
      mesh.jz(i+1, j-1, k)*wtc + 

      mesh.jz(i-1, j  , k)*wts + 
      mesh.jz(i  , j  , k)*wtm + 
      mesh.jz(i+1, j  , k)*wts + 

      mesh.jz(i-1, j+1, k)*wtc + 
      mesh.jz(i  , j+1, k)*wts + 
      mesh.jz(i+1, j+1, k)*wtc;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}


/// 2D 3-point compensator filter from Birdsall & Langdon
//
// Filter is (20,-1,-1) with normalization 1/12
//
//  NOTE: Main difference to a binomial compensator is 
//  the suppressed corners
//
template<>
void fields::Compensator2<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  double winv=1./12.; //normalization
  double wtm=20.0*winv, //middle M
         wts=-1.0*winv, //side   S
         wtc=-1.0*winv; //corner K

  auto& mesh = tile.get_yee();


  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jx(i-1, j-1, k)*wtc + 
      mesh.jx(i  , j-1, k)*wts + 
      mesh.jx(i+1, j-1, k)*wtc + 

      mesh.jx(i-1, j  , k)*wts + 
      mesh.jx(i  , j  , k)*wtm + 
      mesh.jx(i+1, j  , k)*wts + 

      mesh.jx(i-1, j+1, k)*wtc + 
      mesh.jx(i  , j+1, k)*wts + 
      mesh.jx(i+1, j+1, k)*wtc;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jy(i-1, j-1, k)*wtc + 
      mesh.jy(i  , j-1, k)*wts + 
      mesh.jy(i+1, j-1, k)*wtc + 

      mesh.jy(i-1, j  , k)*wts + 
      mesh.jy(i  , j  , k)*wtm + 
      mesh.jy(i+1, j  , k)*wts + 

      mesh.jy(i-1, j+1, k)*wtc + 
      mesh.jy(i  , j+1, k)*wts + 
      mesh.jy(i+1, j+1, k)*wtc;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays

  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
    tmp(i,j,k) = 
      mesh.jz(i-1, j-1, k)*wtc + 
      mesh.jz(i  , j-1, k)*wts + 
      mesh.jz(i+1, j-1, k)*wtc + 

      mesh.jz(i-1, j  , k)*wts + 
      mesh.jz(i  , j  , k)*wtm + 
      mesh.jz(i+1, j  , k)*wts + 

      mesh.jz(i-1, j+1, k)*wtc + 
      mesh.jz(i  , j+1, k)*wts + 
      mesh.jz(i+1, j+1, k)*wtc;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays

}


/// single 2D 3-point general filter pass
template<>
void fields::General3pStrided<2>::solve(
    fields::Tile<2>& tile)
{
  // 2D general coefficients
  double winv=1./4.;                         //normalization
  double wtm=winv * 4.0*alpha*alpha,         //middle
         wts=winv * 2.0*alpha*(1.0-alpha),   //side
         wtc=winv * (1.0-alpha)*(1.0-alpha); //corner

  auto& mesh = tile.get_yee();

  int k = 0;

  //--------------------------------------------------
  // Jx
  //for(int jstr=1; jstr<=stride; jstr++) 
  for(int istr=1; istr<=stride; istr++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jx(i-istr, j-istr, k)*wtc + 
        mesh.jx(i     , j-istr, k)*wts + 
        mesh.jx(i+istr, j-istr, k)*wtc + 

        mesh.jx(i-istr, j     , k)*wts + 
        mesh.jx(i     , j     , k)*wtm + 
        mesh.jx(i+istr, j     , k)*wts + 

        mesh.jx(i-istr, j+istr, k)*wtc + 
        mesh.jx(i     , j+istr, k)*wts + 
        mesh.jx(i+istr, j+istr, k)*wtc;
    }
    mesh.jx = tmp; // then copy from scratch to original arrays
  }

  // Jy
  for(int istr=1; istr<=stride; istr++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jy(i-istr, j-istr, k)*wtc + 
        mesh.jy(i     , j-istr, k)*wts + 
        mesh.jy(i+istr, j-istr, k)*wtc + 
               
        mesh.jy(i-istr, j     , k)*wts + 
        mesh.jy(i     , j     , k)*wtm + 
        mesh.jy(i+istr, j     , k)*wts + 
               
        mesh.jy(i-istr, j+istr, k)*wtc + 
        mesh.jy(i     , j+istr, k)*wts + 
        mesh.jy(i+istr, j+istr, k)*wtc;
    }
    mesh.jy = tmp; // then copy from scratch to original arrays
  }

  // Jz
  for(int istr=1; istr<=stride; istr++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
    for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jz(i-istr, j-istr, k)*wtc + 
        mesh.jz(i     , j-istr, k)*wts + 
        mesh.jz(i+istr, j-istr, k)*wtc + 
               
        mesh.jz(i-istr, j     , k)*wts + 
        mesh.jz(i     , j     , k)*wtm + 
        mesh.jz(i+istr, j     , k)*wts + 
               
        mesh.jz(i-istr, j+istr, k)*wtc + 
        mesh.jz(i     , j+istr, k)*wts + 
        mesh.jz(i+istr, j+istr, k)*wtc;
    }
  mesh.jz = tmp; // then copy from scratch to original arrays
  }

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
  const double wn=1./16.0/16.0;                  //normalization

  auto& mesh = tile.get_yee();

  //--------------------------------------------------
  // Jx
  int k = 0;
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jx(i-3, j-3, 0)*1.0*wn +
        mesh.jx(i-2, j-3, 0)*2.0*wn +
        mesh.jx(i-1, j-3, 0)*3.0*wn +
        mesh.jx(i+0, j-3, 0)*4.0*wn +
        mesh.jx(i+1, j-3, 0)*3.0*wn +
        mesh.jx(i+2, j-3, 0)*2.0*wn +
        mesh.jx(i+3, j-3, 0)*1.0*wn +
        mesh.jx(i-3, j-2, 0)*2.0*wn +
        mesh.jx(i-2, j-2, 0)*4.0*wn +
        mesh.jx(i-1, j-2, 0)*6.0*wn +
        mesh.jx(i+0, j-2, 0)*8.0*wn +
        mesh.jx(i+1, j-2, 0)*6.0*wn +
        mesh.jx(i+2, j-2, 0)*4.0*wn +
        mesh.jx(i+3, j-2, 0)*2.0*wn +
        mesh.jx(i-3, j-1, 0)*3.0*wn +
        mesh.jx(i-2, j-1, 0)*6.0*wn +
        mesh.jx(i-1, j-1, 0)*9.0*wn +
        mesh.jx(i+0, j-1, 0)*12.*wn +
        mesh.jx(i+1, j-1, 0)*9.0*wn +
        mesh.jx(i+2, j-1, 0)*6.0*wn +
        mesh.jx(i+3, j-1, 0)*3.0*wn +
        mesh.jx(i-3, j+0, 0)*4.0*wn +
        mesh.jx(i-2, j+0, 0)*8.0*wn +
        mesh.jx(i-1, j+0, 0)*12.*wn +
        mesh.jx(i+0, j+0, 0)*16.*wn +
        mesh.jx(i+1, j+0, 0)*12.*wn +
        mesh.jx(i+2, j+0, 0)*8.0*wn +
        mesh.jx(i+3, j+0, 0)*4.0*wn +
        mesh.jx(i-3, j+1, 0)*3.0*wn +
        mesh.jx(i-2, j+1, 0)*6.0*wn +
        mesh.jx(i-1, j+1, 0)*9.0*wn +
        mesh.jx(i+0, j+1, 0)*12.*wn +
        mesh.jx(i+1, j+1, 0)*9.0*wn +
        mesh.jx(i+2, j+1, 0)*6.0*wn +
        mesh.jx(i+3, j+1, 0)*3.0*wn +
        mesh.jx(i-3, j+2, 0)*2.0*wn +
        mesh.jx(i-2, j+2, 0)*4.0*wn +
        mesh.jx(i-1, j+2, 0)*6.0*wn +
        mesh.jx(i+0, j+2, 0)*8.0*wn +
        mesh.jx(i+1, j+2, 0)*6.0*wn +
        mesh.jx(i+2, j+2, 0)*4.0*wn +
        mesh.jx(i+3, j+2, 0)*2.0*wn +
        mesh.jx(i-3, j+3, 0)*1.0*wn +
        mesh.jx(i-2, j+3, 0)*2.0*wn +
        mesh.jx(i-1, j+3, 0)*3.0*wn +
        mesh.jx(i+0, j+3, 0)*4.0*wn +
        mesh.jx(i+1, j+3, 0)*3.0*wn +
        mesh.jx(i+2, j+3, 0)*2.0*wn +
        mesh.jx(i+3, j+3, 0)*1.0*wn;
  }
  mesh.jx = tmp; // then copy from scratch to original arrays

  // Jy
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jy(i-3, j-3, 0)*1.0*wn +
        mesh.jy(i-2, j-3, 0)*2.0*wn +
        mesh.jy(i-1, j-3, 0)*3.0*wn +
        mesh.jy(i+0, j-3, 0)*4.0*wn +
        mesh.jy(i+1, j-3, 0)*3.0*wn +
        mesh.jy(i+2, j-3, 0)*2.0*wn +
        mesh.jy(i+3, j-3, 0)*1.0*wn +
        mesh.jy(i-3, j-2, 0)*2.0*wn +
        mesh.jy(i-2, j-2, 0)*4.0*wn +
        mesh.jy(i-1, j-2, 0)*6.0*wn +
        mesh.jy(i+0, j-2, 0)*8.0*wn +
        mesh.jy(i+1, j-2, 0)*6.0*wn +
        mesh.jy(i+2, j-2, 0)*4.0*wn +
        mesh.jy(i+3, j-2, 0)*2.0*wn +
        mesh.jy(i-3, j-1, 0)*3.0*wn +
        mesh.jy(i-2, j-1, 0)*6.0*wn +
        mesh.jy(i-1, j-1, 0)*9.0*wn +
        mesh.jy(i+0, j-1, 0)*12.*wn +
        mesh.jy(i+1, j-1, 0)*9.0*wn +
        mesh.jy(i+2, j-1, 0)*6.0*wn +
        mesh.jy(i+3, j-1, 0)*3.0*wn +
        mesh.jy(i-3, j+0, 0)*4.0*wn +
        mesh.jy(i-2, j+0, 0)*8.0*wn +
        mesh.jy(i-1, j+0, 0)*12.*wn +
        mesh.jy(i+0, j+0, 0)*16.*wn +
        mesh.jy(i+1, j+0, 0)*12.*wn +
        mesh.jy(i+2, j+0, 0)*8.0*wn +
        mesh.jy(i+3, j+0, 0)*4.0*wn +
        mesh.jy(i-3, j+1, 0)*3.0*wn +
        mesh.jy(i-2, j+1, 0)*6.0*wn +
        mesh.jy(i-1, j+1, 0)*9.0*wn +
        mesh.jy(i+0, j+1, 0)*12.*wn +
        mesh.jy(i+1, j+1, 0)*9.0*wn +
        mesh.jy(i+2, j+1, 0)*6.0*wn +
        mesh.jy(i+3, j+1, 0)*3.0*wn +
        mesh.jy(i-3, j+2, 0)*2.0*wn +
        mesh.jy(i-2, j+2, 0)*4.0*wn +
        mesh.jy(i-1, j+2, 0)*6.0*wn +
        mesh.jy(i+0, j+2, 0)*8.0*wn +
        mesh.jy(i+1, j+2, 0)*6.0*wn +
        mesh.jy(i+2, j+2, 0)*4.0*wn +
        mesh.jy(i+3, j+2, 0)*2.0*wn +
        mesh.jy(i-3, j+3, 0)*1.0*wn +
        mesh.jy(i-2, j+3, 0)*2.0*wn +
        mesh.jy(i-1, j+3, 0)*3.0*wn +
        mesh.jy(i+0, j+3, 0)*4.0*wn +
        mesh.jy(i+1, j+3, 0)*3.0*wn +
        mesh.jy(i+2, j+3, 0)*2.0*wn +
        mesh.jy(i+3, j+3, 0)*1.0*wn;
  }
  mesh.jy = tmp; // then copy from scratch to original arrays


  // Jz
  for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) 
  for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
      tmp(i,j,k) = 
        mesh.jz(i-3, j-3, 0)*1.0*wn +
        mesh.jz(i-2, j-3, 0)*2.0*wn +
        mesh.jz(i-1, j-3, 0)*3.0*wn +
        mesh.jz(i+0, j-3, 0)*4.0*wn +
        mesh.jz(i+1, j-3, 0)*3.0*wn +
        mesh.jz(i+2, j-3, 0)*2.0*wn +
        mesh.jz(i+3, j-3, 0)*1.0*wn +
        mesh.jz(i-3, j-2, 0)*2.0*wn +
        mesh.jz(i-2, j-2, 0)*4.0*wn +
        mesh.jz(i-1, j-2, 0)*6.0*wn +
        mesh.jz(i+0, j-2, 0)*8.0*wn +
        mesh.jz(i+1, j-2, 0)*6.0*wn +
        mesh.jz(i+2, j-2, 0)*4.0*wn +
        mesh.jz(i+3, j-2, 0)*2.0*wn +
        mesh.jz(i-3, j-1, 0)*3.0*wn +
        mesh.jz(i-2, j-1, 0)*6.0*wn +
        mesh.jz(i-1, j-1, 0)*9.0*wn +
        mesh.jz(i+0, j-1, 0)*12.*wn +
        mesh.jz(i+1, j-1, 0)*9.0*wn +
        mesh.jz(i+2, j-1, 0)*6.0*wn +
        mesh.jz(i+3, j-1, 0)*3.0*wn +
        mesh.jz(i-3, j+0, 0)*4.0*wn +
        mesh.jz(i-2, j+0, 0)*8.0*wn +
        mesh.jz(i-1, j+0, 0)*12.*wn +
        mesh.jz(i+0, j+0, 0)*16.*wn +
        mesh.jz(i+1, j+0, 0)*12.*wn +
        mesh.jz(i+2, j+0, 0)*8.0*wn +
        mesh.jz(i+3, j+0, 0)*4.0*wn +
        mesh.jz(i-3, j+1, 0)*3.0*wn +
        mesh.jz(i-2, j+1, 0)*6.0*wn +
        mesh.jz(i-1, j+1, 0)*9.0*wn +
        mesh.jz(i+0, j+1, 0)*12.*wn +
        mesh.jz(i+1, j+1, 0)*9.0*wn +
        mesh.jz(i+2, j+1, 0)*6.0*wn +
        mesh.jz(i+3, j+1, 0)*3.0*wn +
        mesh.jz(i-3, j+2, 0)*2.0*wn +
        mesh.jz(i-2, j+2, 0)*4.0*wn +
        mesh.jz(i-1, j+2, 0)*6.0*wn +
        mesh.jz(i+0, j+2, 0)*8.0*wn +
        mesh.jz(i+1, j+2, 0)*6.0*wn +
        mesh.jz(i+2, j+2, 0)*4.0*wn +
        mesh.jz(i+3, j+2, 0)*2.0*wn +
        mesh.jz(i-3, j+3, 0)*1.0*wn +
        mesh.jz(i-2, j+3, 0)*2.0*wn +
        mesh.jz(i-1, j+3, 0)*3.0*wn +
        mesh.jz(i+0, j+3, 0)*4.0*wn +
        mesh.jz(i+1, j+3, 0)*3.0*wn +
        mesh.jz(i+2, j+3, 0)*2.0*wn +
        mesh.jz(i+3, j+3, 0)*1.0*wn;
  }
  mesh.jz = tmp; // then copy from scratch to original arrays


}



//template class fields::Binomial2<1>; // 1D
template class fields::Binomial2<2>; // 2D
template class fields::Binomial2<3>; // 3D


//template class fields::General3p<1>; // 1D
template class fields::General3p<2>; // 2D
//template class fields::General3p<3>; // 3D


//template class fields::General3pStrided<1>; // 1D
template class fields::General3pStrided<2>; // 2D
//template class fields::General3pStrided<3>; // 3D
  
//template class fields::Binomial2Strided2<1>; // 1D
  template class fields::Binomial2Strided2<2>; // 2D
//template class fields::Binomial2Strided2<3>; // 3D
  

//template class fields::Compensator2<1>; // 1D
  template class fields::Compensator2<2>; // 2D
//template class fields::Compensator2<3>; // 3D
