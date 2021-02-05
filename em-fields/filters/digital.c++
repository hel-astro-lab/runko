#include "digital.h"

#include <cmath>


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

  real_short C;

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
    C = C2[is+1][js+1];

    for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,0) += mesh.jx(i+is, j+js, 0)*C;
        }
      }
  }
  swap(mesh.jx, tmp);

  // Jy
  tmp.clear();
  for(int is=-1; is<=1; is++) 
  for(int js=-1; js<=1; js++) {
    C = C2[is+1][js+1];

    for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,0) += mesh.jy(i+is, j+js, 0)*C;
        }
      }
  }
  swap(mesh.jy, tmp);

  // Jz
  tmp.clear();
  for(int is=-1; is<=1; is++) 
  for(int js=-1; js<=1; js++) {
    C = C2[is+1][js+1];

    for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,0) += mesh.jz(i+is, j+js, 0)*C;
        }
      }
  }
  swap(mesh.jz, tmp);
  
}


void sweep_in_x(
        toolbox::Mesh<real_short,3>& arr,
        int lenx,
        int leny,
        int lenz,
        int sub_cycle
        ) {

  int imin, imax, jmin, jmax, kmin, kmax;
  real_short tmp1, tmp2;
    
  for(int npass=1; npass<=sub_cycle; npass++){
    imin = -3 + npass;
    imax = lenx + 3 - 1 - npass;

    jmin = -3 + npass;
    jmax = leny + 3 - 1 - npass;

    kmin = -3 + npass;
    kmax = lenz + 3 - 1 - npass;

    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {

        tmp2 = arr(imin-1, j, k);
        tmp2 = arr(imin-1, j, k);
        for(int i=imin; i<imax-1; i=i+2) {
        //for i = imin, imax - 1, 2 //mod 2 imax-imin
        //for i = imin, imax, 2     //non mod2 

          tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k);

          arr(i - 1, j, k) = tmp2;
          tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i + 1, j, k) + 0.25 * arr(i + 2, j, k);
          arr(i, j, k) = tmp1;
          }

        // extra manual step for uneven grids
        tmp1 = 0.25 * arr(imax - 1, j, k) + 0.5 * arr(imax, j, k) + 0.25 * arr(imax + 1, j, k);
        arr(imax - 1, j, k) = tmp2;

        arr(imax, j, k) = tmp1;
      }//j
    }//k
  } // end of npass
}



/// Single 3D 2nd order 3-point Binomial filter
// NOTE: optimized to slide in x/y/z dimensions to enable vectorization
// TODO: not fully implemented but frame is there
template<>
void fields::OptBinomial2<3>::solve(
    fields::Tile<3>& tile)
{
  auto& mesh = tile.get_yee();

  sweep_in_x(mesh.jx, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);
  sweep_in_x(mesh.jy, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);
  sweep_in_x(mesh.jz, tile.mesh_lengths[0], tile.mesh_lengths[1], tile.mesh_lengths[2], 1);

  //sweep_in_y(mesh.jx, 1);
  //sweep_in_y(mesh.jy, 1);
  //sweep_in_y(mesh.jz, 1);

  //sweep_in_z(mesh.jx, 1);
  //sweep_in_z(mesh.jy, 1);
  //sweep_in_z(mesh.jz, 1);

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

  auto& mesh = tile.get_yee();

  const int halo = 2; // FIXME

  const int imin = 0 - halo;
  const int imax = tile.mesh_lengths[0] + halo;

  const int jmin = 0 - halo;
  const int jmax = tile.mesh_lengths[1] + halo;

  const int kmin = 0 - halo;
  const int kmax = tile.mesh_lengths[2] + halo;


  // NOTE: using tmp as scratch arrays
    
  //--------------------------------------------------
  // Jx
  //NOTE: C3 array varies fastest with nested i(j(k( ...))) looping 

  tmp.clear();
  for(int is=-1; is<=1; is++) {
  for(int js=-1; js<=1; js++) {
  for(int ks=-1; ks<=1; ks++) {
    C = C3[is+1][js+1][ks+1];

    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,k) += mesh.jx(i+is, j+js, k+ks)*C;
        }
      }}
  }}}
  swap(mesh.jx, tmp);

  //--------------------------------------------------
  // Jy
  tmp.clear();
  for(int is=-1; is<=1; is++) {
  for(int js=-1; js<=1; js++) {
  for(int ks=-1; ks<=1; ks++) {
    C = C3[is+1][js+1][ks+1];

    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,k) += mesh.jy(i+is, j+js, k+ks)*C;
        }
      }}
  }}}
  swap(mesh.jy, tmp);

  //--------------------------------------------------
  // Jz
  tmp.clear();
  for(int is=-1; is<=1; is++) {
  for(int js=-1; js<=1; js++) {
  for(int ks=-1; ks<=1; ks++) {
    C = C3[is+1][js+1][ks+1];

    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        #pragma omp simd
        for(int i=imin; i<imax; i++) {
          tmp(i,j,k) += mesh.jz(i+is, j+js, k+ks)*C;
        }
      }}
  }}}
  swap(mesh.jz, tmp);

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
//  NOTE: Main difference to a binomial compensator are
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


/// single 2D 3-point strided filter pass
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
  const double wn=1./16.0/16.0;  //normalization

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
