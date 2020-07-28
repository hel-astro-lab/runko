#include <iostream>
#include <cmath>

#include "tile.h"
#include <nvtx3/nvToolsExt.h> 

namespace fields {


  template<typename F, typename... Args>
  __global__ void
    iterate3D(F fun, int xMax, int yMax, int zMax, Args... args)
  {
    //
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int idz = blockIdx.z * blockDim.z + threadIdx.z;
  
    if(idx >= xMax) return;
    if(idy >= yMax) return;
    if(idz >= zMax) return;
  
    fun(idx, idy, idz, args...);
  }
  

  /*
__global__
void copy_vert_yeeKernel(YeeLattice *lhs, YeeLattice *rhs, int Ny, int Nz, int halo, int ito, int ifro, int in)
{

  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  int h = blockIdx.z * blockDim.z + threadIdx.z;

  if(j >= Ny) return;
  if(k >= Nz) return;
  if(h >= halo) return;

  lhs->ex(ito+in*h, j, k) = rhs->ex(ifro+in*h, j, k);
  lhs->ey(ito+in*h, j, k) = rhs->ey(ifro+in*h, j, k);
  lhs->ez(ito+in*h, j, k) = rhs->ez(ifro+in*h, j, k);

  lhs->bx(ito+in*h, j, k) = rhs->bx(ifro+in*h, j, k);
  lhs->by(ito+in*h, j, k) = rhs->by(ifro+in*h, j, k);
  lhs->bz(ito+in*h, j, k) = rhs->bz(ifro+in*h, j, k);

  lhs->jx(ito+in*h, j, k) = rhs->jx(ifro+in*h, j, k);
  lhs->jy(ito+in*h, j, k) = rhs->jy(ifro+in*h, j, k);
  lhs->jz(ito+in*h, j, k) = rhs->jz(ifro+in*h, j, k);
}
*/

void copy_vert_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Ny, int Nz, int halo, int ito, int ifro, int in)
{
    //

    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 8, 8, 1 };
    dim3 grid  = { 1 + (Ny / 8),
                  1 + (Nz / 8),
                  halo };

    //copy_vert_yeeKernel<<<grid, block>>>(lhs_dev, rhs_dev, Ny, Nz, halo, ito, ifro, in);

    iterate3D<<<grid, block>>>([=] __device__ __host__ (int j, int k ,int h, YeeLattice* lhs_in, YeeLattice* rhs_in){
      //
      lhs_in->ex(ito+in*h, j, k) = rhs_in->ex(ifro+in*h, j, k);
      lhs_in->ey(ito+in*h, j, k) = rhs_in->ey(ifro+in*h, j, k);
      lhs_in->ez(ito+in*h, j, k) = rhs_in->ez(ifro+in*h, j, k);
    
      lhs_in->bx(ito+in*h, j, k) = rhs_in->bx(ifro+in*h, j, k);
      lhs_in->by(ito+in*h, j, k) = rhs_in->by(ifro+in*h, j, k);
      lhs_in->bz(ito+in*h, j, k) = rhs_in->bz(ifro+in*h, j, k);
    
      lhs_in->jx(ito+in*h, j, k) = rhs_in->jx(ifro+in*h, j, k);
      lhs_in->jy(ito+in*h, j, k) = rhs_in->jy(ifro+in*h, j, k);
      lhs_in->jz(ito+in*h, j, k) = rhs_in->jz(ifro+in*h, j, k);
    }, Ny, Nz, halo, lhs_dev, rhs_dev);

    cudaDeviceSynchronize();

}

}