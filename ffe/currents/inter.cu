

#include "rffe.h"
#include "../../tools/signum.h"
#include "../../em-fields/tile.h"

#include <cmath>

#include <iostream>

#include "../../tools/mesh.h"

template <typename F, typename ... Args>
__global__ 
void interateXYZKernVari(F fun, int xMax, int yMax, int zMax, Args ... args)
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

  fun(idx, idy, idz, args...);
}

void interpolateDevEntry( 
  toolbox::Mesh<real_short,3>& f,
  toolbox::Mesh<real_short,0>& fi,
  const std::array<int,3>& in,
  const std::array<int,3>& out
)
{

int im = in[2] == out[2] ? 0 :  -out[2];
int ip = in[2] == out[2] ? 0 : 1-out[2];

int jm = in[1] == out[1] ? 0 :  -out[1];
int jp = in[1] == out[1] ? 0 : 1-out[1];

int km = in[0] == out[0] ? 0 :  -out[0];
int kp = in[0] == out[0] ? 0 : 1-out[0];

dim3 block = {4,4,4};
dim3 grid = {1+(f.Nx/4),1+(f.Ny/4),1+(f.Nz/4)};


cudaHostRegister(&f, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
cudaHostRegister(&fi, sizeof(toolbox::Mesh<real_short,0>), cudaHostRegisterMapped);

toolbox::Mesh<real_short,3> *f_dev;
toolbox::Mesh<real_short,0> *fi_dev;

cudaHostGetDevicePointer(&f_dev, &f, 0);
cudaHostGetDevicePointer(&fi_dev, &fi, 0);

interateXYZKernVari<<<grid, block>>>([=] __device__ __host__ (int i, int j, int k, toolbox::Mesh<real_short,3> *f_in, toolbox::Mesh<real_short,0> *fi_in){
  //
  real_short f11, f10, f01, f00, f1, f0;

  f11 = (*f_in)(i+ip, j+jp, k+km) + (*f_in)(i+ip, j+jp, k+kp);
  f10 = (*f_in)(i+ip, j+jm, k+km) + (*f_in)(i+ip, j+jm, k+kp);
  f01 = (*f_in)(i+im, j+jp, k+km) + (*f_in)(i+im, j+jp, k+kp);
  f00 = (*f_in)(i+im, j+jm, k+km) + (*f_in)(i+im, j+jm, k+kp);
  f1  = f11 + f10;
  f0  = f01 + f00;

  (*fi_in)(i,j,k) = 0.125f*(f1 + f0);
}, f.Nx, f.Ny, f.Nz, f_dev, fi_dev);

/*
auto err = cudaDeviceSynchronize();
if (err != cudaSuccess) 
{
   fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(err));
}
*/
}



void push_ebDevEntry(ffe::Tile<3>& tile)
{
  //
    // refs to storages
    fields::YeeLattice&     m = tile.get_yee();
    ffe::SkinnyYeeLattice& dm = tile.dF; 
  
    // refs to fields for easier access
    /*
    auto& ex  = m.ex;
    auto& ey  = m.ey;
    auto& ez  = m.ez;
  
    auto& bx  = m.bx;
    auto& by  = m.by;
    auto& bz  = m.bz;
  */
    real_short c = tile.cfl;
  
    // dt / dx
    real_short cx = c;
    real_short cy = c;
    real_short cz = c;

    dim3 block = {4,4,4};
    dim3 grid = {1+(static_cast<int>(tile.mesh_lengths[2])/4),1+(static_cast<int>(tile.mesh_lengths[1])/4),1+(static_cast<int>(tile.mesh_lengths[0])/4)};

    cudaHostRegister(&dm, sizeof(ffe::SkinnyYeeLattice), cudaHostRegisterMapped);
    cudaHostRegister(&m, sizeof(fields::YeeLattice), cudaHostRegisterMapped);

    ffe::SkinnyYeeLattice *dm_dev;
    fields::YeeLattice *m_dev;
    cudaHostGetDevicePointer(&dm_dev, &dm, 0);
    cudaHostGetDevicePointer(&m_dev, &m, 0);

    interateXYZKernVari<<<grid, block>>>([=] __device__ __host__ (int i, int j, int k, ffe::SkinnyYeeLattice *dm_in, fields::YeeLattice *m_in){
      //
                // dB = dt*curl E
                dm_in->bx(i,j,k) = cz*( m_in->ey(i,  j,  k+1) - m_in->ey(i,j,k) ) - cy*( m_in->ez(i,  j+1,k) - m_in->ez(i,j,k) );
                dm_in->by(i,j,k) = cx*( m_in->ez(i+1,j,  k  ) - m_in->ez(i,j,k) ) - cz*( m_in->ex(i,  j,k+1) - m_in->ex(i,j,k) );
                dm_in->bz(i,j,k) = cy*( m_in->ex(i,  j+1,k  ) - m_in->ex(i,j,k) ) - cx*( m_in->ey(i+1,j,k  ) - m_in->ey(i,j,k) );
        
                // dE = dt*curl B 
                dm_in->ex(i,j,k) = cz*( m_in->by(i,  j,  k-1) - m_in->by(i,j,k) ) - cy*( m_in->bz( i,  j-1,k  ) - m_in->bz(i,j,k) );
                dm_in->ey(i,j,k) = cx*( m_in->bz(i-1,j,  k  ) - m_in->bz(i,j,k) ) - cz*( m_in->bx( i,  j,  k-1) - m_in->bx(i,j,k) );
                dm_in->ez(i,j,k) = cy*( m_in->bx(i,  j-1,k  ) - m_in->bx(i,j,k) ) - cx*( m_in->by( i-1,j,  k  ) - m_in->by(i,j,k) );
    }, static_cast<int>(tile.mesh_lengths[2]),
    static_cast<int>(tile.mesh_lengths[1]),
    static_cast<int>(tile.mesh_lengths[0]), dm_dev, m_dev);
 
    /*
    auto err = cudaDeviceSynchronize();
if (err != cudaSuccess) 
{
   fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(err));
}
*/
}


template<>
void ffe::rFFE2<3>::add_jperpXDevEntry(ffe::Tile<3>& tile)
{
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  auto& jx  = m.jx;
  //auto& jy  = m.jy;
  //auto& jz  = m.jz;

  real_short dt = tile.cfl;

  dim3 block = {4,4,4};
  dim3 grid = {1+(static_cast<int>(tile.mesh_lengths[2])/4),1+(static_cast<int>(tile.mesh_lengths[1])/4),1+(static_cast<int>(tile.mesh_lengths[0])/4)};

  cudaHostRegister(&dm, sizeof(ffe::SkinnyYeeLattice), cudaHostRegisterMapped);
  cudaHostRegister(&jx, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  cudaHostRegister(&bxf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&byf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&bzf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&exf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&eyf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&ezf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&rhf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  ffe::SkinnyYeeLattice *dm_dev;
  toolbox::Mesh<real_short,3> *jx_dev;

  toolbox::Mesh<real_short, 0> *bxf_dev;
  toolbox::Mesh<real_short, 0> *byf_dev;
  toolbox::Mesh<real_short, 0> *bzf_dev;
  toolbox::Mesh<real_short, 0> *exf_dev;
  toolbox::Mesh<real_short, 0> *eyf_dev;
  toolbox::Mesh<real_short, 0> *ezf_dev;
  toolbox::Mesh<real_short, 0> *rhf_dev;


  cudaHostGetDevicePointer(&dm_dev, &dm, 0);
  cudaHostGetDevicePointer(&jx_dev, &jx, 0);

  cudaHostGetDevicePointer(&bxf_dev, &bxf, 0);
  cudaHostGetDevicePointer(&byf_dev, &byf, 0);
  cudaHostGetDevicePointer(&bzf_dev, &bzf, 0);
  cudaHostGetDevicePointer(&exf_dev, &exf, 0);
  cudaHostGetDevicePointer(&eyf_dev, &eyf, 0);
  cudaHostGetDevicePointer(&ezf_dev, &ezf, 0);
  cudaHostGetDevicePointer(&rhf_dev, &rhf, 0);

  //https://developer.nvidia.com/blog/new-compiler-features-cuda-8/
  interateXYZKernVari<<<grid, block>>>([=] __device__ __host__ (int i, int j, int k, ffe::SkinnyYeeLattice *dm_in, toolbox::Mesh<real_short,3> *jx_in,
  toolbox::Mesh<real_short, 0> *bxf_in,
  toolbox::Mesh<real_short, 0> *byf_in,
  toolbox::Mesh<real_short, 0> *bzf_in,
  toolbox::Mesh<real_short, 0> *exf_in,
  toolbox::Mesh<real_short, 0> *eyf_in,
  toolbox::Mesh<real_short, 0> *ezf_in,
  toolbox::Mesh<real_short, 0> *rhf_in
  ){
  real_short b2, cur;
        b2 = (
            (*bxf_in)(i,j,k)*(*bxf_in)(i,j,k) 
          + (*byf_in)(i,j,k)*(*byf_in)(i,j,k) 
          + (*bzf_in)(i,j,k)*(*bzf_in)(i,j,k) 
          + EPS);

        cur = (*rhf_in)(i,j,k) * ((*eyf_in)(i,j,k)*(*bzf_in)(i,j,k) - (*byf_in)(i,j,k)*(*ezf_in)(i,j,k) )/b2;
        (*jx_in)(i,j,k) = cur;
        dm_in->ex(i,j,k) -= dt*cur;
  }, static_cast<int>(tile.mesh_lengths[2]),
  static_cast<int>(tile.mesh_lengths[1]),
  static_cast<int>(tile.mesh_lengths[0]), dm_dev, jx_dev, bxf_dev,
  byf_dev,
  bzf_dev,
  exf_dev,
  eyf_dev,
  ezf_dev,
  rhf_dev);
/*
  auto err = cudaDeviceSynchronize();
  if (err != cudaSuccess) 
  {
     fprintf(stderr,"GPUassert: %s add_jperpXDevEntry %d\n", cudaGetErrorString(err));
  }
  */
}

template<>
void ffe::rFFE2<3>::add_jperpYDevEntry(ffe::Tile<3>& tile)
{
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  //auto& jx  = m.jx;
  auto& jy  = m.jy;
  //auto& jz  = m.jz;

  real_short dt = tile.cfl;

  dim3 block = {4,4,4};
  dim3 grid = {1+(static_cast<int>(tile.mesh_lengths[2])/4),1+(static_cast<int>(tile.mesh_lengths[1])/4),1+(static_cast<int>(tile.mesh_lengths[0])/4)};

  cudaHostRegister(&dm, sizeof(ffe::SkinnyYeeLattice), cudaHostRegisterMapped);
  cudaHostRegister(&jy, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  cudaHostRegister(&bxf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&byf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&bzf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&exf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&eyf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&ezf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&rhf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  ffe::SkinnyYeeLattice *dm_dev;
  toolbox::Mesh<real_short,3> *jy_dev;

  toolbox::Mesh<real_short, 0> *bxf_dev;
  toolbox::Mesh<real_short, 0> *byf_dev;
  toolbox::Mesh<real_short, 0> *bzf_dev;
  toolbox::Mesh<real_short, 0> *exf_dev;
  toolbox::Mesh<real_short, 0> *eyf_dev;
  toolbox::Mesh<real_short, 0> *ezf_dev;
  toolbox::Mesh<real_short, 0> *rhf_dev;

  cudaHostGetDevicePointer(&dm_dev, &dm, 0);
  cudaHostGetDevicePointer(&jy_dev, &jy, 0);

  cudaHostGetDevicePointer(&bxf_dev, &bxf, 0);
  cudaHostGetDevicePointer(&byf_dev, &byf, 0);
  cudaHostGetDevicePointer(&bzf_dev, &bzf, 0);
  cudaHostGetDevicePointer(&exf_dev, &exf, 0);
  cudaHostGetDevicePointer(&eyf_dev, &eyf, 0);
  cudaHostGetDevicePointer(&ezf_dev, &ezf, 0);
  cudaHostGetDevicePointer(&rhf_dev, &rhf, 0);

  //https://developer.nvidia.com/blog/new-compiler-features-cuda-8/
  interateXYZKernVari<<<grid, block>>>([=] __device__ __host__ (int i, int j, int k, ffe::SkinnyYeeLattice *dm_in, toolbox::Mesh<real_short,3> *jy_in,
  toolbox::Mesh<real_short, 0> *bxf_in,
  toolbox::Mesh<real_short, 0> *byf_in,
  toolbox::Mesh<real_short, 0> *bzf_in,
  toolbox::Mesh<real_short, 0> *exf_in,
  toolbox::Mesh<real_short, 0> *eyf_in,
  toolbox::Mesh<real_short, 0> *ezf_in,
  toolbox::Mesh<real_short, 0> *rhf_in
  ){
  real_short b2, cur;
  b2 = (
    (*bxf_in)(i,j,k)*(*bxf_in)(i,j,k) 
    + (*byf_in)(i,j,k)*(*byf_in)(i,j,k) 
    + (*bzf_in)(i,j,k)*(*bzf_in)(i,j,k) 
    + EPS);

cur = (*rhf_in)(i,j,k) * ((*ezf_in)(i,j,k)*(*bxf_in)(i,j,k) - (*exf_in)(i,j,k)*(*bzf_in)(i,j,k))/b2;
        (*jy_in)(i,j,k) = cur;
        dm_in->ey(i,j,k) -= dt*cur;
  }, static_cast<int>(tile.mesh_lengths[2]),
  static_cast<int>(tile.mesh_lengths[1]),
  static_cast<int>(tile.mesh_lengths[0]), dm_dev, jy_dev, bxf_dev,
  byf_dev,
  bzf_dev,
  exf_dev,
  eyf_dev,
  ezf_dev,
  rhf_dev);
  
}

template<>
void ffe::rFFE2<3>::add_jperpZDevEntry(ffe::Tile<3>& tile)
{
  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  //auto& jx  = m.jx;
  //auto& jy  = m.jy;
  auto& jz  = m.jz;

  real_short dt = tile.cfl;

  dim3 block = {4,4,4};
  dim3 grid = {1+(static_cast<int>(tile.mesh_lengths[2])/4),1+(static_cast<int>(tile.mesh_lengths[1])/4),1+(static_cast<int>(tile.mesh_lengths[0])/4)};

  cudaHostRegister(&dm, sizeof(ffe::SkinnyYeeLattice), cudaHostRegisterMapped);
  cudaHostRegister(&jz, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  cudaHostRegister(&bxf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&byf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&bzf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&exf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&eyf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&ezf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);
  cudaHostRegister(&rhf, sizeof(toolbox::Mesh<real_short,3>), cudaHostRegisterMapped);

  ffe::SkinnyYeeLattice *dm_dev;
  toolbox::Mesh<real_short,3> *jz_dev;

  toolbox::Mesh<real_short, 0> *bxf_dev;
  toolbox::Mesh<real_short, 0> *byf_dev;
  toolbox::Mesh<real_short, 0> *bzf_dev;
  toolbox::Mesh<real_short, 0> *exf_dev;
  toolbox::Mesh<real_short, 0> *eyf_dev;
  toolbox::Mesh<real_short, 0> *ezf_dev;
  toolbox::Mesh<real_short, 0> *rhf_dev;

  cudaHostGetDevicePointer(&dm_dev, &dm, 0);
  cudaHostGetDevicePointer(&jz_dev, &jz, 0);

  cudaHostGetDevicePointer(&bxf_dev, &bxf, 0);
  cudaHostGetDevicePointer(&byf_dev, &byf, 0);
  cudaHostGetDevicePointer(&bzf_dev, &bzf, 0);
  cudaHostGetDevicePointer(&exf_dev, &exf, 0);
  cudaHostGetDevicePointer(&eyf_dev, &eyf, 0);
  cudaHostGetDevicePointer(&ezf_dev, &ezf, 0);
  cudaHostGetDevicePointer(&rhf_dev, &rhf, 0);

  //https://developer.nvidia.com/blog/new-compiler-features-cuda-8/
  interateXYZKernVari<<<grid, block>>>([=] __device__ __host__ (int i, int j, int k, ffe::SkinnyYeeLattice *dm_in, toolbox::Mesh<real_short,3> *jz_in,
  toolbox::Mesh<real_short, 0> *bxf_in,
  toolbox::Mesh<real_short, 0> *byf_in,
  toolbox::Mesh<real_short, 0> *bzf_in,
  toolbox::Mesh<real_short, 0> *exf_in,
  toolbox::Mesh<real_short, 0> *eyf_in,
  toolbox::Mesh<real_short, 0> *ezf_in,
  toolbox::Mesh<real_short, 0> *rhf_in
  ){
  real_short b2, cur;
  b2 = (
    (*bxf_in)(i,j,k)*(*bxf_in)(i,j,k) 
    + (*byf_in)(i,j,k)*(*byf_in)(i,j,k) 
    + (*bzf_in)(i,j,k)*(*bzf_in)(i,j,k) 
    + EPS);

cur = (*rhf_in)(i,j,k) * ((*exf_in)(i,j,k)*(*byf_in)(i,j,k) - (*bxf_in)(i,j,k)*(*eyf_in)(i,j,k))/b2;
        (*jz_in)(i,j,k) = cur;
        dm_in->ez(i,j,k) -= dt*cur;
  }, static_cast<int>(tile.mesh_lengths[2]),
  static_cast<int>(tile.mesh_lengths[1]),
  static_cast<int>(tile.mesh_lengths[0]), dm_dev, jz_dev, bxf_dev,
  byf_dev,
  bzf_dev,
  exf_dev,
  eyf_dev,
  ezf_dev,
  rhf_dev);
  
}




template<>
void ffe::rFFE2<3>::remove_jparDevEntry(ffe::Tile<3>& tile)
{
  //

  nvtxRangePush(__FUNCTION__);

  fields::YeeLattice&     m = tile.get_yee();
  ffe::SkinnyYeeLattice& dm = tile.dF; 

  real_short cur, b2;
  real_short dt = tile.cfl;

  stagger_x_eb(m);
  cudaDeviceSynchronize();
  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bxf(i,j,k) /b2/dt;

        m.jx(i,j,k) += cur;
        //dm.ex(i,j,k) = m.ex(i,j,k) - cur;
          
        //m.ex(i,j,k) -= cur*dt;
        dm.ex(i,j,k) = m.ex(i,j,k) - cur*dt;
      }
    }
  }

  stagger_y_eb(m);
  cudaDeviceSynchronize();
  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*byf(i,j,k) /b2/dt;

        m.jy(i,j,k) += cur;
        //dm.ey(i,j,k) = m.ey(i,j,k) - cur;
        //m.ey(i,j,k) -= cur*dt;
        
        dm.ey(i,j,k) = m.ey(i,j,k) - cur*dt;
      }
    }
  }

  stagger_z_eb(m);
  cudaDeviceSynchronize();
  for(int k=0; k<static_cast<int>(tile.mesh_lengths[2]); k++) {
    for(int j=0; j<static_cast<int>(tile.mesh_lengths[1]); j++) {
      for(int i=0; i<static_cast<int>(tile.mesh_lengths[0]); i++) {
        b2 = (
            bxf(i,j,k)*bxf(i,j,k) + 
            byf(i,j,k)*byf(i,j,k) + 
            bzf(i,j,k)*bzf(i,j,k) + 
            EPS);
        cur = (exf(i,j,k)*bxf(i,j,k) + eyf(i,j,k)*byf(i,j,k) + ezf(i,j,k)*bzf(i,j,k))*bzf(i,j,k) /b2/dt;

        m.jz(i,j,k) += cur;
        //dm.ez(i,j,k) = m.ez(i,j,k) - cur;
        //m.ez(i,j,k) -= cur*dt;

        dm.ez(i,j,k) = m.ez(i,j,k) - cur*dt;
      }
    }
  }

nvtxRangePop();
}