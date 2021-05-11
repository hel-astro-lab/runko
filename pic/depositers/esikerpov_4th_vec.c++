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

  const int vec_size = 16; 
  unsigned int batch_size = 7*7*7*vec_size; // = 22 Mbytes (16 vec)

  float batch_J[ batch_size ] __attribute__( ( aligned( 64 ) ) );
  float Sx0_v[48]             __attribute__( ( aligned( 64 ) ) );
  float Sy0_v[48]             __attribute__( ( aligned( 64 ) ) );
  float Sz0_v[48]             __attribute__( ( aligned( 64 ) ) );
  float DSx_v[56]             __attribute__( ( aligned( 64 ) ) );
  float DSy_v[56]             __attribute__( ( aligned( 64 ) ) );
  float DSz_v[56]             __attribute__( ( aligned( 64 ) ) );


  static constexpr float one_third    = 1.0/3.0;
  static constexpr float f_1_ov_384   = 1.0/384.0;
  static constexpr float f_1_ov_48    = 1.0/48.0;
  static constexpr float f_1_ov_16    = 1.0/16.0;
  static constexpr float f_1_ov_12    = 1.0/12.0;
  static constexpr float f_1_ov_24    = 1.0/24.0;
  static constexpr float f_19_ov_96   = 19.0/96.0;
  static constexpr float f_11_ov_24   = 11.0/24.0;
  static constexpr float f_1_ov_4     = 1.0/4.0;
  static constexpr float f_1_ov_6     = 1.0/6.0;
  static constexpr float f_115_ov_192 = 115.0/192.0;
  static constexpr float f_5_ov_8     = 5.0/8.0;


  for(auto&& con : tile.containers) {

    int istart = 0;
    int iend   = con.size();
    int nparts = iend-istart;


    // Closest multiple of vec_size higher or equal to npart = iend-istart.
    //int num_vecs = ( iend-istart+(nparts-1)-( (iend-istart-1)&(nparts-1)))/vec_size;
    //if(num_vecs*vec_size != nparts) num_vecs++;

    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge


    #pragma omp simd
    for( unsigned int j=0; j<batch_size; j++ ) batch_J[j] = 0.;


    // loop over batches of particles; each batch holds vec_size x particles
    // that are processed simultaneously.
    for(int ivect=0; ivect < nparts; ivect += vec_size) {

      // process either vec_size prtcls or up to the end of array
      int np_computed = min(nparts-ivect, vec_size);

      #pragma omp simd
      for( int ip=0 ; ip<np_computed; ip++ ) {

        //--------------------------------------------------
        // particle locations
        
        int n = ip + ivect*nparts; // actual running prtcl loop index
        float u = con.vel(0,n);
        float v = con.vel(1,n);
        float w = con.vel(2,n);
        float invgam = 1.0/sqrt(1.0 + u*u + v*v + w*w);
          
        // new (normalized) location, x_{n+1}
        float x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
        float y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
        float z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

        // previous location, x_n
        float x1 = x2 - u*invgam*c;
        float y1 = y2 - v*invgam*c;
        float z1 = z2 - w*invgam*c; 

        // primary grid; -1/2 to +1/2
        int i1p = D >= 1 ? round(x1) : 0;
        int i2p = D >= 1 ? round(x2) : 0;
        int j1p = D >= 2 ? round(y1) : 0;
        int j2p = D >= 2 ? round(y2) : 0;
        int k1p = D >= 3 ? round(z1) : 0;
        int k2p = D >= 3 ? round(z2) : 0;

        //--------------------------------------------------
        // calculate shape coefficients at former time step
        float d;
          
        d = x1 - float(i1p);
        Sx0_v[            ip] = f1_ov_384   - f1_ov_48  * d  + f1_ov_16 * d*d - f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sx0_v[1*vec_size +ip] = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sx0_v[2*vec_size +ip] = f115_ov_192                  - f5_ov_8  * d*d                    + f1_ov_4  * d*d*d*d;
        Sx0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sx0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sx0_v[5*vec_size +ip] = 0.0f;
        
        //                             Y                                 //
        d = y1 - float(j1p);
        Sy0_v[            ip] = f1_ov_384   - f1_ov_48  * d  + f1_ov_16 * d*d - f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sy0_v[1*vec_size +ip] = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sy0_v[2*vec_size +ip] = f115_ov_192                  - f5_ov_8  * d*d                    + f1_ov_4  * d*d*d*d;
        Sy0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sy0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sy0_v[5*vec_size +ip] = 0.0f;
        
        //                             Z                                 //
        d = z1 - float(k1p);
        Sz0_v[            ip] = f1_ov_384   - f1_ov_48  * d  + f1_ov_16 * d*d - f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sz0_v[1*vec_size +ip] = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sz0_v[2*vec_size +ip] = f115_ov_192                  - f5_ov_8  * d*d                    + f1_ov_4  * d*d*d*d;
        Sz0_v[3*vec_size +ip] = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        Sz0_v[4*vec_size +ip] = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        Sz0_v[5*vec_size +ip] = 0.0f;

        //--------------------------------------------------
        // calculate shape coefficients at new time step

        float S0, S1, S2, S3, S4, m1, c0, p1; // temporary variables
        int cell_shift;

        d = x2 - float(i2p), 
        cell_shift = i2p - i1p; 

        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S1 = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S2 = f115_ov_192                  - f5_ov_8  * d*d                    + f1_ov_4  * d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;

        DSx_v[            ip] = m1*S0;
        DSx_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sx0_v[           ip];
        DSx_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sx0_v[  vec_size+ip];
        DSx_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sx0_v[2*vec_size+ip];
        DSx_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sx0_v[3*vec_size+ip];
        DSx_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sx0_v[4*vec_size+ip];
        DSx_v[6*vec_size +ip] =                                 p1*S4;


        //--------------------------------------------------
        // y
        d = y2 - float(j2p), 
        cell_shift = j2p - j1p; 
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S0 = f1_ov_384   - f1_ov_48  * d  + f1_ov_16 * d*d - f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        S1 = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S2 = f115_ov_192                  - f5_ov_8  * d*d                    + f1_ov_4  * d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;

        DSy_v[            ip] = m1*S0;
        DSy_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sy0_v[           ip];
        DSy_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sy0_v[  vec_size+ip];
        DSy_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sy0_v[2*vec_size+ip];
        DSy_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sy0_v[3*vec_size+ip];
        DSy_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sy0_v[4*vec_size+ip];
        DSy_v[6*vec_size +ip] =                                 p1*S4;


        //--------------------------------------------------
        // z
        d = z2 - float(k2p), 
        cell_shift = k2p - k1p; 
        m1 = ( cell_shift == -1 );
        c0 = ( cell_shift ==  0 );
        p1 = ( cell_shift ==  1 );

        S0 = f1_ov_384   - f1_ov_48  * d  + f1_ov_16 * d*d - f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;
        S1 = f19_ov_96   - f11_ov_24 * d  + f1_ov_4  * d*d + f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S2 = f115_ov_192 - f5_ov_8   * d*d + f1_ov_4  * d*d*d*d;
        S3 = f19_ov_96   + f11_ov_24 * d  + f1_ov_4  * d*d - f1_ov_6  * d*d*d - f1_ov_6  * d*d*d*d;
        S4 = f1_ov_384   + f1_ov_48  * d  + f1_ov_16 * d*d + f1_ov_12 * d*d*d + f1_ov_24 * d*d*d*d;

        DSz_v[            ip] = m1*S0;
        DSz_v[  vec_size +ip] = c0*S0 + m1*S1                         - Sz0_v[           ip];
        DSz_v[2*vec_size +ip] = p1*S0 + c0*S1 + m1*S2                 - Sz0_v[  vec_size+ip];
        DSz_v[3*vec_size +ip] =         p1*S1 + c0*S2 + m1*S3         - Sz0_v[2*vec_size+ip];
        DSz_v[4*vec_size +ip] =                 p1*S2 + c0*S3 + m1*S4 - Sz0_v[3*vec_size+ip];
        DSz_v[5*vec_size +ip] =                         p1*S3 + c0*S4 - Sz0_v[4*vec_size+ip];
        DSz_v[6*vec_size +ip] =                                    p1 * S4;
    

      //} // end if ip // FIXME not in GPU
      //// Jx^(d,p,p)
      //#pragma omp simd
      //for( int ip=0 ; ip<np_computed; ip++ ) {

        //--------------------------------------------------
        //// Jx^(d,p,p)
        float sum[7] = {0.0f};
        for(int k=1; k<7; k++) sum[k] = sum[k-1]-DSx_v[( k-1 )*vec_size+ip];
        
        float tmp = q*c*one_third*DSy_v[ip]*DSz_v[ip];
        for(int i=1; i<7; i++) batch_J [ i*49*vec_size + ip ] += sum[i]*tmp;
        

        for(int k=1 ; k<7 ; k++ ) {
          tmp = q*c*(0.5f*DSy_v[ip]*Sz0_buff_vect[ (k-1)*vec_size + ip] + one_third*DSy_v[ip]*DSz_v[k*vec_size + ip] );
          int index = k*vec_size + ip;

          for(int i=1; i<7; i++) batch_J[ index + 49*i*vec_size ] += sum[i]*tmp;
        }

        for(int j=1 ; j<7 ; j++ ) {
          tmp = crx_p * ( 0.5f*DSz_v[ip]*Sy0_buff_vect[( j-1 )*vec_size+ip] + one_third*DSy_v[j*vec_size+ip]*DSz_v[ip] );
          int index = j*7*vec_size+ip;

          for(int i=1; i<7; i++) batch_J[ index+49*i*vec_size ] += sum[i]*tmp;
        }

        for( int j=1 ; j<7 ; j++ ) {
          for( int k=1 ; k<7 ; k++ ) {
            tmp = q*c*( Sy0_buff_vect[( j-1 )*vec_size+ip]*Sz0_buff_vect[( k-1 )*vec_size+ip]
                                + 0.5f*DSy_v[j*vec_size+ip]*Sz0_buff_vect[( k-1 )*vec_size+ip]
                                + 0.5f*DSz_v[k*vec_size+ip]*Sy0_buff_vect[( j-1 )*vec_size+ip]
                                + one_third*DSy_v[j*vec_size+ip]*DSz_v[k*vec_size+ip] );
            int index = (j*7 + k)*vec_size + ip;
            for(int i=1; i<7; i++) batch_J[ index + 49*i*vec_size ] += sum[i]*tmp;
          }
        }


    // FIXME not in GPU
    //  } // end of ip
    //} // end of ivect

    //-------------------------------------------------- 
    // deposit cached jx current to grid
    int iloc0 = ipom2*nprimy*nprimz+jpom2*nprimz+kpom2;

    int iloc  = iloc0;
    for( int i=1 ; i<7 ; i++ ) {
      iloc += nprimy*nprimz;
      for( int j=0 ; j<7 ; j++ ) {

        #pragma omp simd
        for( int k=0 ; k<7 ; k++ ) {
          float tmpJx = 0.;
          int ilocal = ( ( i )*49+j*7+k )*vec_size;

          #pragma unroll(8)
          for( int ip=0 ; ip<8; ip++ ) tmpJx += batch_J[ilocal + ip];

          //Jx[iloc+j*nprimz+k] += tmpJx;
          //atomic_add( yee.jx(iloc, jloc, kloc), tmpJx[j][k] );
          yee.jx(iloc, jloc, kloc) += tmpJx[j][k];
        }
      }
    }

    // TODO: why all calcs are repeated???
    //  why not just compute these inside the loop?

    //-------------------------------------------------- 
    //-------------------------------------------------- 
    //-------------------------------------------------- 
    // Jy^(p,d,p)

    #pragma omp simd
    for( int j=0; j<batch_size; j++ ) batch_J[j] = 0.0f;



    nparts = ( int )iend-( int )istart;
    for( int ivect=0 ; ivect < nparts; ivect += vec_size ) {
      int np_computed( min( nparts-ivect, vec_size ) );
        

      #pragma omp simd
      for( int ip=0 ; ip<np_computed; ip++ ) {
        float sum[7] = {0.0f};

        for( int k=1 ; k<7 ; k++ ) sum[k] = sum[k-1]-DSy[( k-1 )*vec_size+ip];
        
        float tmp( cry_p *one_third*DSz[ip]*DSx[ip] );
        for( int j=1 ; j<7 ; j++ ) batch_J [( ( j )*7 )*vec_size+ip] += sum[j]*tmp;


        for( int k=1 ; k<7 ; k++ ) {
            tmp = cry_p * ( 0.5*DSx[0]*Sz0_buff_vect[( k-1 )*vec_size+ip] + one_third*DSz[k*vec_size+ip]*DSx[ip] );
            int index( ( k )*vec_size+ip );

            for( int j=1 ; j<7 ; j++ ) batch_J [ index+7*j*vec_size ] += sum[j]*tmp;
        }
        for( int i=1 ; i<7 ; i++ ) {
            tmp = cry_p * ( 0.5*DSz[ip]*Sx0_buff_vect[( i-1 )*vec_size+ip] + one_third*DSz[0]*DSx[i*vec_size+ip] );
            int index( ( i*49 )*vec_size+ip );

            for( int j=1 ; j<7 ; j++ ) batch_J [ index+7*j*vec_size ] += sum[j]*tmp;
        }
        for( int i=1 ; i<7 ; i++ ) {
            for( int k=1 ; k<7 ; k++ ) {
                tmp = cry_p * ( Sz0_buff_vect[( k-1 )*vec_size+ip]*Sx0_buff_vect[( i-1 )*vec_size+ip]
                                + 0.5*DSz[k*vec_size+ip]*Sx0_buff_vect[( i-1 )*vec_size+ip]
                                + 0.5*DSx[i*vec_size+ip]*Sz0_buff_vect[( k-1 )*vec_size+ip]
                                + one_third*DSz[k*vec_size+ip]*DSx[i*vec_size+ip] );
                int index( ( i*49 + k )*vec_size+ip );

                for( int j=1 ; j<7 ; j++ ) batch_J [ index+7*j*vec_size ] += sum[j]*tmp;
            }
        }//i
      }
    } // end of ivect


    iloc = iloc0+ipom2*nprimz;
    for( int i=0 ; i<7 ; i++ ) {
        for( int j=1 ; j<7 ; j++ ) {
            #pragma omp simd
            for( int k=0 ; k<7 ; k++ ) {
                float tmpJy = 0.;
                int ilocal = ( ( i )*49+j*7+k )*vec_size;

                #pragma unroll(8)
                for( int ip=0 ; ip<8; ip++ ) tmpJy += batch_J [ilocal+ip];

                Jy[iloc+j*nprimz+k] += tmpJy;
            }
        }
        iloc += ( nprimy+1 )*nprimz;
    }


    //-------------------------------------------------- 
    //-------------------------------------------------- 
    //-------------------------------------------------- 
    // Jz^(p,p,d)

    nparts = ( int )iend-( int )istart;
    #pragma omp simd
    for( int j=0; j<batch_size; j++ ) batch_J[j] = 0.;


    for( int ivect=0 ; ivect < nparts; ivect += vec_size ) {
        int np_computed( min( nparts-ivect, vec_size ) );


        #pragma omp simd
        for( int ip=0 ; ip<np_computed; ip++ ) {
            
            float sum[7] = {0.0f};
            for( int k=1 ; k<7 ; k++ ) sum[k] = sum[k-1]-DSz[( k-1 )*vec_size+ip];
            

            float tmp = q*c*one_third*DSx[ip]*DSy[ip];
            for( int k=1 ; k<7 ; k++ ) batch_J[( k )*vec_size+ip] += sum[k]*tmp;

            for( int j=1 ; j<7 ; j++ ) {
                tmp = q*c*( 0.5f*DSx[ip]*Sy0_buff_vect[( j-1 )*vec_size+ip] + one_third*DSx[ip]*DSy[j*vec_size+ip] );

                int index = j*7*vec_size + ip;
                for( int k=1 ; k<7 ; k++ ) batch_J [ index+k*vec_size ] += sum[k]*tmp;
            }

            for( int i=1 ; i<7 ; i++ ) {
                tmp = q*c*( 0.5f*DSy[ip]*Sx0_buff_vect[( i-1 )*vec_size+ip] + one_third*DSx[i*vec_size+ip]*DSy[ip] );
                int index = (i*49)*vec_size + ip;

                for( int k=1 ; k<7 ; k++ ) batch_J [ index+k*vec_size ] += sum[k]*tmp;
            }

            for( int i=1 ; i<7 ; i++ ) {
                for( int j=1 ; j<7 ; j++ ) {
                    tmp = q*c*d( Sx0_buff_vect[( i-1 )*vec_size+ip]*Sy0_buff_vect[( j-1 )*vec_size+ip]
                                    + 0.5f*DSx[i*vec_size+ip]*Sy0_buff_vect[( j-1 )*vec_size+ip]
                                    + 0.5f*DSy[j*vec_size+ip]*Sx0_buff_vect[( i-1 )*vec_size+ip]
                                    + one_third*DSx[i*vec_size+ip]*DSy[j*vec_size+ip] );
                    int index( ( i*49 + j*7 )*vec_size+ip );

                    for( int k=1 ; k<7 ; k++ ) batch_J [ index+k*vec_size ] += sum[k]*tmp;
                }
            }
            
        } // end of ip
    }// end of ivect



    iloc = iloc0  + jpom2 +ipom2*nprimy;
    for( int i=0 ; i<7 ; i++ ) {
        for( int j=0 ; j<7 ; j++ ) {

            #pragma omp simd
            for( int k=1 ; k<7 ; k++ ) {
                float tmpJz = 0.;
                int ilocal = ( ( i )*49+j*7+k )*vec_size;

                #pragma unroll(8)
                for( int ip=0 ; ip<8; ip++ ) tmpJz +=  batch_J[ilocal+ip];

                Jz [iloc + ( j )*( nprimz+1 ) + k] +=  tmpJz;
            }
        }
        iloc += nprimy*( nprimz+1 );
    }



  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
//template class pic::Esikerpov_2nd<1,3>; // 1D3V
//template class pic::Esikerpov_2nd<2,3>; // 2D3V
template class pic::Esikerpov_4th<3,3>; // 3D3V

