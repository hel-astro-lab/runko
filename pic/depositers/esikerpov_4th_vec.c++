#include "esikerpov_4th.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "../shapes.h"
#include "../../tools/iter/iter.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;

template<typename T>
inline T clamp(T x, T xmin, T xmax) {
  x = std::max(xmin, x);
  x = std::min(xmax, x);
  return x;
}

// capped index so that it does not overflow tile from max side; 
// stencil offset is baked in via the -3 terms; N+2 is maximum of tile with halo
size_t cap(int ip, size_t i, int N) {
  //return i + size_t(clamp(ip+int(i)-3, -3, N+2)) - size_t(ip); 
  return 3 + std::min(ip+int(i)-3, N+2) - ip;
}



template<size_t D, size_t V>
void pic::Esikerpov_4th<D,V>::solve( pic::Tile<D>& tile )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& yee = tile.get_yee();
  const auto mins = tile.mins;

  //clear arrays before new update
  yee.jx.clear();
  yee.jy.clear();
  yee.jz.clear();
  yee.rho.clear();

  const auto Nx = yee.Nx;
  const auto Ny = yee.Ny;
  const auto Nz = yee.Nz;

  //--------------------------------------------------
  // vectorized arrays

  const int vec_size = 8; 
  //unsigned int batch_size = 7*7*7*vec_size; // = 22 Mbytes (16 vec)
  unsigned int bs = 7*vec_size;

  //double batch_J[ batch_size ] __attribute__( ( aligned( 64 ) ) );
  double Sx0_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double Sy0_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double Sz0_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double DSx_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double DSy_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double DSz_v[bs]  __attribute__( ( aligned( 64 ) ) );
  double sum_v[bs]  __attribute__( ( aligned( 64 ) ) );


  static constexpr double one_third   = 1.0/3.0;
  static constexpr double f1_ov_384   = 1.0/384.0;
  static constexpr double f1_ov_48    = 1.0/48.0;
  static constexpr double f1_ov_16    = 1.0/16.0;
  static constexpr double f1_ov_12    = 1.0/12.0;
  static constexpr double f1_ov_24    = 1.0/24.0;
  static constexpr double f19_ov_96   = 19.0/96.0;
  static constexpr double f11_ov_24   = 11.0/24.0;
  static constexpr double f1_ov_4     = 1.0/4.0;
  static constexpr double f1_ov_6     = 1.0/6.0;
  static constexpr double f115_ov_192 = 115.0/192.0;
  static constexpr double f5_ov_8     = 5.0/8.0;


  for(auto&& con : tile.containers) {

    int istart = 0;
    int iend   = con.size();
    int nparts = iend-istart;

    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? yee.jx.indx(0,1,0) - yee.jx.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? yee.jx.indx(0,0,1) - yee.jx.indx(0,0,0) : 0;

    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge


    // loop over batches of particles; each batch holds vec_size x particles
    // that are processed simultaneously.
    //
    // NOTE: ivect is the index of the first prtcl of the batch
    for(int ivect=0; ivect < nparts; ivect += vec_size) {

      // process either vec_size prtcls or up to the end of array
      int batch_size = min(nparts-ivect, vec_size);

      #pragma omp simd
      for(int ip=0 ; ip<batch_size; ip++ ) {

        //--------------------------------------------------
        // particle locations
        
        int n = ip + ivect; // actual running prtcl loop index
        double u = con.vel(0,n);
        double v = con.vel(1,n);
        double w = con.vel(2,n);
        double invgam = 1.0/sqrt(1.0 + u*u + v*v + w*w);
          
        // new (normalized) location, x_{n+1}
        double x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
        double y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
        double z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

        // previous location, x_n
        double x1 = x2 - u*invgam*c;
        double y1 = y2 - v*invgam*c;
        double z1 = z2 - w*invgam*c; 

        // primary grid; -1/2 to +1/2
        int i1p = D >= 1 ? round(x1) : 0;
        int j1p = D >= 2 ? round(y1) : 0;
        int k1p = D >= 3 ? round(z1) : 0;

        int i2p = D >= 1 ? round(x2) : 0;
        int j2p = D >= 2 ? round(y2) : 0;
        int k2p = D >= 3 ? round(z2) : 0;


        //--------------------------------------------------
        // calculate shape coefficients at former time step
        double d;
          
        d = x1 - double(i1p);
        Sx0_v[            ip] = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sx0_v[  vec_size +ip] = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sx0_v[2*vec_size +ip] = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        Sx0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sx0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sx0_v[5*vec_size +ip] = 0.0;

        d = y1 - double(j1p);
        Sy0_v[            ip] = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sy0_v[  vec_size +ip] = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sy0_v[2*vec_size +ip] = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        Sy0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sy0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sy0_v[5*vec_size +ip] = 0.0;
        
        d = z1 - double(k1p);
        Sz0_v[            ip] = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sz0_v[  vec_size +ip] = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sz0_v[2*vec_size +ip] = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        Sz0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        Sz0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        Sz0_v[5*vec_size +ip] = 0.0;


        //--------------------------------------------------
        // calculate shape coefficients at new time step

        double S0, S1, S2, S3, S4, m1, c0, p1; // temporary variables
        int cell_shift;

        d = x2 - double(i2p), 
        cell_shift = i2p - i1p; 

        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S0 = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        S1 = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S2 = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;

        DSx_v[            ip] = m1*S0;
        DSx_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sx0_v[           ip];
        DSx_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sx0_v[  vec_size+ip];
        DSx_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sx0_v[2*vec_size+ip];
        DSx_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sx0_v[3*vec_size+ip];
        DSx_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sx0_v[4*vec_size+ip];
        DSx_v[6*vec_size +ip] =                                 p1*S4;


        //--------------------------------------------------
        // y
        d = y2 - double(j2p), 
        cell_shift = j2p - j1p; 
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S0 = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        S1 = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S2 = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;

        DSy_v[            ip] = m1*S0;
        DSy_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sy0_v[           ip];
        DSy_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sy0_v[  vec_size+ip];
        DSy_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sy0_v[2*vec_size+ip];
        DSy_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sy0_v[3*vec_size+ip];
        DSy_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sy0_v[4*vec_size+ip];
        DSy_v[6*vec_size +ip] =                                 p1*S4;


        //--------------------------------------------------
        // z
        d = z2 - double(k2p), 
        cell_shift = k2p - k1p; 
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S0 = f1_ov_384   - f1_ov_48 *d  + f1_ov_16*d*d - f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;
        S1 = f19_ov_96   - f11_ov_24*d  + f1_ov_4 *d*d + f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S2 = f115_ov_192                - f5_ov_8 *d*d                  + f1_ov_4 *d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24*d  + f1_ov_4 *d*d - f1_ov_6 *d*d*d - f1_ov_6 *d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48 *d  + f1_ov_16*d*d + f1_ov_12*d*d*d + f1_ov_24*d*d*d*d;

        DSz_v[            ip] = m1*S0;
        DSz_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sz0_v[           ip];
        DSz_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sz0_v[  vec_size+ip];
        DSz_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sz0_v[2*vec_size+ip];
        DSz_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sz0_v[3*vec_size+ip];
        DSz_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sz0_v[4*vec_size+ip];
        DSz_v[6*vec_size +ip] =                                 p1*S4;
    

        //std::cout << "index1" << i1p << " " << j1p << " " << k1p 
        //          << "index2" << i2p << " " << j2p << " " << k2p << "\n";
        //std::cout 
        //  << "cur " 
        //  << " i1: "  << "(" << i1p << "," << j1p << "," << k1p << ")"
        //  << " i2: "  << "(" << i2p << "," << j2p << "," << k2p << ")"
        //  << " x1: " << "(" << x1 << "," << y1 << "," << z1 << ")"
        //  << " x2: " << "(" << x2 << "," << y2 << "," << z2 << ")"
        //  << "\n";

        //size_t ind = yee.jx.indx(i2p,j2p,k2p);  // TODO wrong
        //std::cout << "done" << ind << "\n";

        //size_t ind = yee.jx.indx(i1p,j1p,k1p); 


        //--------------------------------------------------
        //// Jx^(d,p,p)
        sum_v[ip+0*vec_size] = 0.0;
        // FIXME changed
        for(size_t k=1; k<7; k++) sum_v[ip + k*vec_size] = sum_v[ip+(k-1)*vec_size] - q*c*DSx_v[ip+(k-1)*vec_size];

        //int  z_size0 = nprimz;
        //int yz_size0 = nprimz*nprimy;
        //int linindex0 = iold[ipart+0*nparts]*yz_size0+iold[ipart+1*nparts]*z_size0+iold[ipart+2*nparts];

        // one-dimensional index
        size_t ind0 = yee.jx.indx(i1p,j1p,k1p); 

        for(size_t k=0; k<7; k++ ) {
          size_t k0 = cap(k1p, k, Nz);

          for(size_t j=0; j<7; j++ ) {
            size_t j0 = cap(j1p, j, Ny);

            double tmp =             Sy0_v[ip + j*vec_size]*Sz0_v[ip + k*vec_size] 
                         +       0.5*DSy_v[ip + j*vec_size]*Sz0_v[ip + k*vec_size] 
                         +       0.5*DSz_v[ip + k*vec_size]*Sy0_v[ip + j*vec_size] 
                         + one_third*DSy_v[ip + j*vec_size]*DSz_v[ip + k*vec_size];

            for(size_t i=1; i<7; i++) {
              size_t i0 = cap(i1p, i, Nx);

              double val = sum_v[ip + i*vec_size]*tmp;
              size_t idx = ind0 + (i0-3) + (j0-3)*iy + (k0-3)*iz;

              assert(val < 1.0);

              atomic_add( yee.jx(idx), val );
              //yee.jx(idx) += val;
            }
          }
        }


        //--------------------------------------------------
        // Jy^(p,d,p)
        sum_v[ip+0*vec_size] = 0.0;
        for(size_t k=1; k<7; k++) sum_v[ip+k*vec_size] = sum_v[ip+(k-1)*vec_size] - q*c*DSy_v[ip+(k-1)*vec_size];

        //int  z_size1 = nprimz;
        //int yz_size1 = nprimz*( nprimy+1 ); // TODO: has nprimy+1
        //int linindex1 = iold[ipart+0*nparts]*yz_size1+iold[ipart+1*nparts]*z_size1+iold[ipart+2*nparts];

        size_t ind1 = yee.jx.indx(i1p,j1p+1,k1p); 

        for(size_t k=0; k<7; k++ ) {
          size_t k0 = cap(k1p, k, Nz);

          for(size_t i=0; i<7; i++ ) {
            size_t i0 = cap(i1p, i, Nx);

            double tmp =             Sz0_v[ip+k*vec_size]*Sx0_v[ip+i*vec_size] 
                         +       0.5*DSz_v[ip+k*vec_size]*Sx0_v[ip+i*vec_size] 
                         +       0.5*DSx_v[ip+i*vec_size]*Sz0_v[ip+k*vec_size] 
                         + one_third*DSz_v[ip+k*vec_size]*DSx_v[ip+i*vec_size];

            for(size_t j=1; j<7; j++ ) {
              size_t j0 = cap(j1p, j, Ny);
              //size_t idx = ind + i0 + j0*iy + k0*iz; //linindex1 + k + i*yz_size1;
                
              double val = sum_v[ip + j*vec_size]*tmp;
              //size_t idx = ind + (i0-3) + (j0-3)*iy + (k0-3)*iz;
              size_t idx = ind1 + (i0-3) + (j0-3)*(iy+1) + (k0-3)*iz; // FIXME

              assert(val < 1.0);
              atomic_add( yee.jy(idx), val ); // NOTE: jumps of j*iy
              //yee.jy(idx) += val;
            }
          }
        }


        //--------------------------------------------------
        // Jz^(p,p,d)
        sum_v[ip+0*vec_size] = 0.0;
        for(size_t k=1; k<7; k++) sum_v[ip+k*vec_size] = sum_v[ip+(k-1)*vec_size] - q*c*DSz_v[ip+(k-1)*vec_size];

        //int z_size2 =  nprimz+1; // TODO has nprimz+1
        //int yz_size2 = ( nprimz+1 )*nprimy; // TODO has nprimz+1
        //int linindex2 = iold[ipart+0*nparts]*yz_size2+iold[ipart+1*nparts]*z_size2+iold[ipart+2*nparts];

        size_t ind2 = yee.jx.indx(i1p,j1p,k1p+1); 
        for(size_t k=1; k<7; k++ ) {
          size_t k0 = cap(k1p, k, Nz);

          for(size_t i=0; i<7; i++ ) {
            size_t i0 = cap(i1p, i, Nx);

            for(size_t j=0; j<7; j++ ) {
              size_t j0 = cap(j1p, j, Ny);

              double tmp =             Sx0_v[ip+i*vec_size]*Sy0_v[ip+j*vec_size] 
                           +       0.5*DSx_v[ip+i*vec_size]*Sy0_v[ip+j*vec_size] 
                           +       0.5*DSy_v[ip+j*vec_size]*Sx0_v[ip+i*vec_size] 
                           + one_third*DSx_v[ip+i*vec_size]*DSy_v[ip+j*vec_size];

              //size_t idx = linindex2 + j*z_size2 + i*yz_size2;

              double val = sum_v[ip + k*vec_size]*tmp;
              //size_t idx = ind + (i0-3) + (j0-3)*iy + (k0-3)*(iz);
              size_t idx = ind2 + (i0-3) + (j0-3)*iy + (k0-3)*(iz+1); // FIXME

              //std::cout << sum_v[ip + k*vec_size] << " " << tmp;
              assert(val < 1.0);
              atomic_add( yee.jz(idx), val ); // NOTE: jumps of k*iz
              //yee.jz(idx) += val; 
            }
          }
        }


      }// end ip
    } //end ivect
  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
//template class pic::Esikerpov_2nd<1,3>; // 1D3V
//template class pic::Esikerpov_2nd<2,3>; // 2D3V
template class pic::Esikerpov_4th<3,3>; // 3D3V

