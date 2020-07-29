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


}


void copy_horz_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Nz, int halo, int jto, int jfro, int jn)
{

    //

    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 8, 8, 1 };
    dim3 grid  = { 1 + (Nx / 8),
                  1 + (Nz / 8),
                  halo };

    //copy_vert_yeeKernel<<<grid, block>>>(lhs_dev, rhs_dev, Ny, Nz, halo, ito, ifro, in);

    iterate3D<<<grid, block>>>([=] __device__ __host__ (int i, int k ,int g, YeeLattice* lhs_in, YeeLattice* rhs_in){
      lhs_in->ex(i, jto+jn*g, k) = rhs_in->ex(i, jfro+jn*g, k);
      lhs_in->ey(i, jto+jn*g, k) = rhs_in->ey(i, jfro+jn*g, k);
      lhs_in->ez(i, jto+jn*g, k) = rhs_in->ez(i, jfro+jn*g, k);

      lhs_in->bx(i, jto+jn*g, k) = rhs_in->bx(i, jfro+jn*g, k);
      lhs_in->by(i, jto+jn*g, k) = rhs_in->by(i, jfro+jn*g, k);
      lhs_in->bz(i, jto+jn*g, k) = rhs_in->bz(i, jfro+jn*g, k);

      lhs_in->jx(i, jto+jn*g, k) = rhs_in->jx(i, jfro+jn*g, k);
      lhs_in->jy(i, jto+jn*g, k) = rhs_in->jy(i, jfro+jn*g, k);
      lhs_in->jz(i, jto+jn*g, k) = rhs_in->jz(i, jfro+jn*g, k);
    }, Nx, Nz, halo, lhs_dev, rhs_dev);

}


void copy_face_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int Ny, int halo, int kto, int kfro, int kn)
{
    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 8, 8, 1 };
    dim3 grid  = { 1 + (Nx / 8), 1 + (Ny / 8), halo };

    iterate3D<<<grid, block>>>([=] __device__ __host__ (int i, int j ,int g, YeeLattice* lhs_in, YeeLattice* rhs_in){  
        lhs_in->ex(i, j, kto+kn*g) = rhs_in->ex(i, j, kfro+kn*g);
        lhs_in->ey(i, j, kto+kn*g) = rhs_in->ey(i, j, kfro+kn*g);
        lhs_in->ez(i, j, kto+kn*g) = rhs_in->ez(i, j, kfro+kn*g);

        lhs_in->bx(i, j, kto+kn*g) = rhs_in->bx(i, j, kfro+kn*g);
        lhs_in->by(i, j, kto+kn*g) = rhs_in->by(i, j, kfro+kn*g);
        lhs_in->bz(i, j, kto+kn*g) = rhs_in->bz(i, j, kfro+kn*g);

        lhs_in->jx(i, j, kto+kn*g) = rhs_in->jx(i, j, kfro+kn*g);
        lhs_in->jy(i, j, kto+kn*g) = rhs_in->jy(i, j, kfro+kn*g);
        lhs_in->jz(i, j, kto+kn*g) = rhs_in->jz(i, j, kfro+kn*g);
    }, Nx, Ny, halo, lhs_dev, rhs_dev);

}



void copy_x_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nx, int halo, int jto, int jfro, int kto, int kfro, int jn, int kn)
{
    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 32, 1, 1 };
    dim3 grid  = { 1 + (Nx / 32), halo, halo };

    
    iterate3D<<<grid, block>>>([=] __device__ __host__ (int i, int g ,int h, YeeLattice* lhs_in, YeeLattice* rhs_in){
      lhs_in->ex(i, jto+jn*h, kto+kn*g) = rhs_in->ex(i, jfro+jn*h, kfro+kn*g);
      lhs_in->ey(i, jto+jn*h, kto+kn*g) = rhs_in->ey(i, jfro+jn*h, kfro+kn*g);
      lhs_in->ez(i, jto+jn*h, kto+kn*g) = rhs_in->ez(i, jfro+jn*h, kfro+kn*g);

      lhs_in->bx(i, jto+jn*h, kto+kn*g) = rhs_in->bx(i, jfro+jn*h, kfro+kn*g);
      lhs_in->by(i, jto+jn*h, kto+kn*g) = rhs_in->by(i, jfro+jn*h, kfro+kn*g);
      lhs_in->bz(i, jto+jn*h, kto+kn*g) = rhs_in->bz(i, jfro+jn*h, kfro+kn*g);

      lhs_in->jx(i, jto+jn*h, kto+kn*g) = rhs_in->jx(i, jfro+jn*h, kfro+kn*g);
      lhs_in->jy(i, jto+jn*h, kto+kn*g) = rhs_in->jy(i, jfro+jn*h, kfro+kn*g);
      lhs_in->jz(i, jto+jn*h, kto+kn*g) = rhs_in->jz(i, jfro+jn*h, kfro+kn*g);
    },  Nx, halo, halo, lhs_dev, rhs_dev);
  }

  void copy_y_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Ny, int halo, int ito, int ifro, int kto, int kfro, int in, int kn)
{
    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 32, 1, 1 };
    dim3 grid  = { 1 + (Ny / 32), halo, halo };

    
    iterate3D<<<grid, block>>>([=] __device__ __host__ (int j, int g ,int h, YeeLattice* lhs_in, YeeLattice* rhs_in){

      lhs_in->ex(ito+in*h, j, kto+kn*g) = rhs_in->ex(ifro+in*h, j, kfro+kn*g);
      lhs_in->ey(ito+in*h, j, kto+kn*g) = rhs_in->ey(ifro+in*h, j, kfro+kn*g);
      lhs_in->ez(ito+in*h, j, kto+kn*g) = rhs_in->ez(ifro+in*h, j, kfro+kn*g);

      lhs_in->bx(ito+in*h, j, kto+kn*g) = rhs_in->bx(ifro+in*h, j, kfro+kn*g);
      lhs_in->by(ito+in*h, j, kto+kn*g) = rhs_in->by(ifro+in*h, j, kfro+kn*g);
      lhs_in->bz(ito+in*h, j, kto+kn*g) = rhs_in->bz(ifro+in*h, j, kfro+kn*g);

      lhs_in->jx(ito+in*h, j, kto+kn*g) = rhs_in->jx(ifro+in*h, j, kfro+kn*g);
      lhs_in->jy(ito+in*h, j, kto+kn*g) = rhs_in->jy(ifro+in*h, j, kfro+kn*g);
      lhs_in->jz(ito+in*h, j, kto+kn*g) = rhs_in->jz(ifro+in*h, j, kfro+kn*g);
    },  Ny, halo, halo, lhs_dev, rhs_dev);
  }

  void copy_z_pencil_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int Nz, int halo, int ito, int ifro, int jto, int jfro, int in, int jn)
{
    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { 32, 1, 1 };
    dim3 grid  = { 1 + (Nz / 32), halo, halo };

    
    iterate3D<<<grid, block>>>([=] __device__ __host__ (int k, int g ,int h, YeeLattice* lhs_in, YeeLattice* rhs_in){
      lhs_in->ex(ito+in*h, jto+jn*g, k) = rhs_in->ex(ifro+in*h, jfro+jn*g, k);
      lhs_in->ey(ito+in*h, jto+jn*g, k) = rhs_in->ey(ifro+in*h, jfro+jn*g, k);
      lhs_in->ez(ito+in*h, jto+jn*g, k) = rhs_in->ez(ifro+in*h, jfro+jn*g, k);

      lhs_in->bx(ito+in*h, jto+jn*g, k) = rhs_in->bx(ifro+in*h, jfro+jn*g, k);
      lhs_in->by(ito+in*h, jto+jn*g, k) = rhs_in->by(ifro+in*h, jfro+jn*g, k);
      lhs_in->bz(ito+in*h, jto+jn*g, k) = rhs_in->bz(ifro+in*h, jfro+jn*g, k);

      lhs_in->jx(ito+in*h, jto+jn*g, k) = rhs_in->jx(ifro+in*h, jfro+jn*g, k);
      lhs_in->jy(ito+in*h, jto+jn*g, k) = rhs_in->jy(ifro+in*h, jfro+jn*g, k);
      lhs_in->jz(ito+in*h, jto+jn*g, k) = rhs_in->jz(ifro+in*h, jfro+jn*g, k);
  
    },  Nz, halo, halo, lhs_dev, rhs_dev);
  }



  void copy_point_yeeDevEntry(YeeLattice& lhs, YeeLattice& rhs, int halo, int ito, int ifro, int jto, int jfro, int kto, int kfro, int in, int jn, int kn)
  {
    cudaHostRegister(&lhs, sizeof(YeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&rhs, sizeof(YeeLattice), cudaHostRegisterMapped);
  
    YeeLattice *lhs_dev;
    YeeLattice *rhs_dev;
    cudaHostGetDevicePointer(&lhs_dev, &lhs, 0);
    cudaHostGetDevicePointer(&rhs_dev, &rhs, 0);

    dim3 block = { halo, halo, halo };
    dim3 grid  = { 1,1,1 };


    iterate3D<<<grid, block>>>([=] __device__ __host__ (int f, int g ,int h, YeeLattice* lhs_in, YeeLattice* rhs_in){
      lhs_in->ex(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->ex(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->ey(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->ey(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->ez(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->ez(ifro+in*h, jfro+jn*g, kfro+kn*f);

      lhs_in->bx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->bx(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->by(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->by(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->bz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->bz(ifro+in*h, jfro+jn*g, kfro+kn*f);

      lhs_in->jx(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->jx(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->jy(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->jy(ifro+in*h, jfro+jn*g, kfro+kn*f);
      lhs_in->jz(ito +in*h, jto +jn*g, kto +kn*f) =  rhs_in->jz(ifro+in*h, jfro+jn*g, kfro+kn*f);
    },  halo, halo, halo, lhs_dev, rhs_dev);
  }
}