

#include "rffe.h"
#include "../../tools/signum.h"
#include "../../em-fields/tile.h"

#include <cmath>

#include <iostream>

#include "../../tools/mesh.h"

__host__ __device__
int index(int i, int j, int k, int Nx, int Ny, int H)
{
  int indx = (i + H) + (Nx + 2*H)*( (j + H) + (Ny + 2*H)*(k + H));
  return indx;
}


template <typename F>
__global__ 
void interateXYZKern(int xMax, int yMax, int zMax, toolbox::Mesh<real_short,3> &f, toolbox::Mesh<real_short,0> &fi, F fun)
{
  //
  int idx = blockIdx.x *blockDim.x + threadIdx.x;
  int idy = blockIdx.y *blockDim.y + threadIdx.y;
  int idz = blockIdx.z *blockDim.z + threadIdx.z;

  if(idx >= xMax)
    return;
  if(idy >= yMax)
    return;
  if(idz >= zMax)
    return;

  fun(idx, idy, idz, f, fi);
}

__global__ 
void interpolateDevKern(int im, int ip, int jm, int jp, int km, int kp, int nX, int nY, int nZ, real_short *f, real_short *fi)
{
  //
  int idx = blockIdx.x *blockDim.x + threadIdx.x;
  int idy = blockIdx.y *blockDim.y + threadIdx.y;
  int idz = blockIdx.z *blockDim.z + threadIdx.z;

  if(idx >= nX)
    return;
  if(idy >= nY)
    return;
  if(idz >= nZ)
    return;

  int i = idx;
  int j = idy;
  int k = idz;
  //
  real_short f11, f10, f01, f00, f1, f0;

  f11 = f[index(i+ip, j+jp, k+km, nX, nY, 3)] + f[index(i+ip, j+jp, k+kp, nX, nY, 3)];
  f10 = f[index(i+ip, j+jm, k+km, nX, nY, 3)] + f[index(i+ip, j+jm, k+kp, nX, nY, 3)];
  f01 = f[index(i+im, j+jp, k+km, nX, nY, 3)] + f[index(i+im, j+jp, k+kp, nX, nY, 3)];
  f00 = f[index(i+im, j+jm, k+km, nX, nY, 3)] + f[index(i+im, j+jm, k+kp, nX, nY, 3)];
  f1  = f11 + f10;
  f0  = f01 + f00;

  fi[index(i,j,k, nX, nY, 0)] = 0.125f*(f1 + f0);

}



void interpolateDevEntry( 
  toolbox::Mesh<real_short,3>& f,
  toolbox::Mesh<real_short,0>& fi,
  const std::array<int,3>& in,
  const std::array<int,3>& out
)
{
//std::cout << "calling iterpolate "<< f.Nx << " " << f.Ny << " " << f.Nz << " " << std::endl;
int im = in[2] == out[2] ? 0 :  -out[2];
int ip = in[2] == out[2] ? 0 : 1-out[2];

int jm = in[1] == out[1] ? 0 :  -out[1];
int jp = in[1] == out[1] ? 0 : 1-out[1];

int km = in[0] == out[0] ? 0 :  -out[0];
int kp = in[0] == out[0] ? 0 : 1-out[0];

dim3 block = {4,4,4};
dim3 grid = {1+(f.Nx/4),1+(f.Ny/4),1+(f.Nz/4)};


interateXYZKern<<<grid, block>>>(f.Nx, f.Ny, f.Nz, f, fi,[=] __device__ __host__ (int i, int j, int k, toolbox::Mesh<real_short,3> &f_in, toolbox::Mesh<real_short,0> &fi_in){
  //
  real_short f11, f10, f01, f00, f1, f0;

  f11 = f_in(i+ip, j+jp, k+km) + f_in(i+ip, j+jp, k+kp);
  f10 = f_in(i+ip, j+jm, k+km) + f_in(i+ip, j+jm, k+kp);
  f01 = f_in(i+im, j+jp, k+km) + f_in(i+im, j+jp, k+kp);
  f00 = f_in(i+im, j+jm, k+km) + f_in(i+im, j+jm, k+kp);
  f1  = f11 + f10;
  f0  = f01 + f00;

  fi_in(i,j,k) = 0.125f*(f1 + f0);
});

//interpolateDevKern<<<grid, block>>>(im, ip, jm, jp, km, kp, f.Nx, f.Ny, f.Nz, f.data(), fi.data());
/*
real_short f11, f10, f01, f00, f1, f0;

for(int k=0; k<f.Nz; k++) {
for(int j=0; j<f.Ny; j++) {
for(int i=0; i<f.Nx; i++) {
  f11 = f(i+ip, j+jp, k+km) + f(i+ip, j+jp, k+kp);
  f10 = f(i+ip, j+jm, k+km) + f(i+ip, j+jm, k+kp);
  f01 = f(i+im, j+jp, k+km) + f(i+im, j+jp, k+kp);
  f00 = f(i+im, j+jm, k+km) + f(i+im, j+jm, k+kp);
  f1  = f11 + f10;
  f0  = f01 + f00;

  fi(i,j,k) = 0.125*(f1 + f0);
}
}
}
*/
}