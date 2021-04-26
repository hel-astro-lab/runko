#include "zigzag_2nd.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;


// Triangular 2nd order charge fluxes
inline void W2nd(float_m x, float_m xr, int i, float_m* out)
{
    // W_{i+1)^(1) 
    float_m xm = 0.5f*(x + xr) - i; 
    out[0] = xm;
      
    // charge fluxes for triangular shapes
    out[1] = 0.50f*(0.5f - xm)*(0.5f - xm); //W2_im1 
    out[2] = 0.75f - xm*xm;                 //W2_i   
    out[3] = 0.50f*(0.5f + xm)*(0.5f + xm); //W2_ip1 
}


template<size_t D, size_t V>
void pic::ZigZag_2nd<D,V>::solve( pic::Tile<D>& tile )
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


  for(auto&& con : tile.containers) {

    const double c = tile.cfl;    // speed of light
    const float_m q = con.q; // charge

    //UniIter::iterate([=] DEVCALLABLE (
    //            size_t n, 
    //            fields::YeeLattice &yee,
    //            pic::ParticleContainer<D>& con
    //            ){
    for(size_t n=0; n<con.size(); n++) {

      //--------------------------------------------------
      double u = con.vel(0,n);
      double v = con.vel(1,n);
      double w = con.vel(2,n);
      double invgam = 1.0/sqrt(1.0 + u*u + v*v + w*w);

      //--------------------------------------------------
      // new (normalized) location, x_{n+1}
        
      //float_m x1, y1, z1;
      float_m x2 = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      float_m y2 = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      float_m z2 = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      // previos location, x_n
      float_m x1 = x2 - u*invgam*c;
      float_m y1 = y2 - v*invgam*c;
      float_m z1 = z2 - w*invgam*c; 

      //--------------------------------------------------
      int i1  = D >= 1 ? floor(x1) : 0;
      int i2  = D >= 1 ? floor(x2) : 0;
      int j1  = D >= 2 ? floor(y1) : 0;
      int j2  = D >= 2 ? floor(y2) : 0;
      int k1  = D >= 3 ? floor(z1) : 0;
      int k2  = D >= 3 ? floor(z2) : 0;

      // 1st order relay point; +1 is equal to +\Delta x
      //float_m xr = min( float_m(min(i1,i2)+1), max( float_m(max(i1,i2)), float_m(0.5*(x1+x2)) ) );
      //float_m yr = min( float_m(min(j1,j2)+1), max( float_m(max(j1,j2)), float_m(0.5*(y1+y2)) ) );
      //float_m zr = min( float_m(min(k1,k2)+1), max( float_m(max(k1,k2)), float_m(0.5*(z1+z2)) ) );
        
      // 2nd order relay point; +1 is equal to +\Delta x
      float_m xr = min( float_m(min(i1,i2)+1), max( float_m(i1+i2)*0.5f, float_m(0.5f*(x1+x2)) ) );
      float_m yr = min( float_m(min(j1,j2)+1), max( float_m(j1+j2)*0.5f, float_m(0.5f*(y1+y2)) ) );
      float_m zr = min( float_m(min(k1,k2)+1), max( float_m(k1+k2)*0.5f, float_m(0.5f*(z1+z2)) ) );

      //--------------------------------------------------
      // particle weights 
      float_m Wxx1[4], Wxx2[4], Wyy1[4], Wyy2[4], Wzz1[4], Wzz2[4];

      W2nd(x1, xr, i1, Wxx1);
      W2nd(x2, xr, i2, Wxx2);

      W2nd(y1, yr, j1, Wyy1);
      W2nd(y2, yr, j2, Wyy2);

      W2nd(z1, zr, k1, Wzz1);
      W2nd(z2, zr, k2, Wzz2);

      //reduce dimension if needed
      //if(D <= 1) {
      //  Wy1a, Wy2a, Wy3a = 1.0;
      //  Wy1b, Wy2b, Wy3b = 1.0;
      //}
      //if(D <= 2) {
      //  Wz1a, Wz2a, Wz3a = 1.0;
      //  Wz1b, Wz2b, Wz3b = 1.0;
      //}


      //Wx1_ip1 = 0.5f*(x + xr) - i; 
      //jx(i-1, j-1, k) = qvx * (0.5 - Wx1_ip1)*W2_jm1
      //jx(i  , j-1, k) = qvx * (0.5 + Wx1_ip1)*W2_jm1
      //jx(i-1, j  , k) = qvx * (0.5 - Wx1_ip1)*W2_j  
      //jx(i  , j  , k) = qvx * (0.5 + Wx1_ip1)*W2_j  
      //jx(i-1, j+1, k) = qvx * (0.5 - Wx1_ip1)*W2_jp1
      //jx(i  , j+1, k) = qvx * (0.5 + Wx1_ip1)*W2_jp1

      //--------------------------------------------------
      // q v = q (x_{i+1} - x_i)/dt
      //
      // NOTE: +q since - sign is already included in the Ampere's equation
      // NOTE: extra c to introduce time step immediately; therefore we store on grid J -> J\Delta t
      // NOTE: More generally we should have: q = weight*qe;
      float_m qvx1 = +q*(xr - x1);
      float_m qvy1 = +q*(yr - y1);
      float_m qvz1 = +q*(zr - z1);

      float_m qvx2 = +q*(x2 - xr);
      float_m qvy2 = +q*(y2 - yr);
      float_m qvz2 = +q*(z2 - zr);
        
      //current deposited at xmid = 0.5(x2-x1)
      //or think this is the relay pt
      //in its own frame prtcl prtcl footpoitn is same
      //in lab frame ftpnt is compressed by 1/gam
      //edge at i-1 experiences inv gam less charge
      //edge at i   experiences inv gam less charge


      //--------------------------------------------------
        
      //jx
      for(int zi=-1; zi <=1; zi++)
      for(int yi=-1; yi <=1; yi++){
        //std::cout << "jx: injecting into" <<
        //"(" << i1-1 <<","<< j1+yi <<","<< k1+zi <<") " <<
        //"(" << i1   <<","<< j1+yi <<","<< k1+zi <<") " <<
        //"(" << i2-1 <<","<< j2+yi <<","<< k2+zi <<") " <<
        //"(" << i2   <<","<< j2+yi <<","<< k2+zi <<") " << "\n";

        //first part of trajectory
        atomic_add( yee.jx(i1-1, j1+yi, k1+zi), qvx1*(0.5f - Wxx1[0])*Wyy1[yi+2]*Wzz1[zi+2] );
        atomic_add( yee.jx(i1  , j1+yi, k1+zi), qvx1*(0.5f + Wxx1[0])*Wyy1[yi+2]*Wzz1[zi+2] );

        //second part of trajectory
        atomic_add( yee.jx(i2-1, j2+yi, k2+zi), qvx2*(0.5f - Wxx2[0])*Wyy2[yi+2]*Wzz2[zi+2] );
        atomic_add( yee.jx(i2,   j2+yi, k2+zi), qvx2*(0.5f + Wxx2[0])*Wyy2[yi+2]*Wzz2[zi+2] );
      }

      //jy
      for(int zi=-1; zi <=1; zi++)
      for(int xi=-1; xi <=1; xi++){
        //std::cout << "jy: injecting into" <<
        //"(" << i1+xi <<","<< j1-1 <<","<< k1+zi <<") " <<
        //"(" << i1+xi <<","<< j1   <<","<< k1+zi <<") " <<
        //"(" << i2+xi <<","<< j2-1 <<","<< k2+zi <<") " <<
        //"(" << i2+xi <<","<< j2   <<","<< k2+zi <<") " << "\n";

        atomic_add( yee.jy(i1+xi, j1-1, k1+zi), qvy1*(0.5f - Wyy1[0])*Wxx1[xi+2]*Wzz1[zi+2] );
        atomic_add( yee.jy(i1+xi, j1,   k1+zi), qvy1*(0.5f + Wyy1[0])*Wxx1[xi+2]*Wzz1[zi+2] );

        atomic_add( yee.jy(i2+xi, j2-1, k2+zi), qvy2*(0.5f - Wyy2[0])*Wxx2[xi+2]*Wzz2[zi+2] );
        atomic_add( yee.jy(i2+xi, j2,   k2+zi), qvy2*(0.5f + Wyy2[0])*Wxx2[xi+2]*Wzz2[zi+2] );
      }

      //jz
      for(int yi=-1; yi <=1; yi++)
      for(int xi=-1; xi <=1; xi++){
        //std::cout << "jz: injecting into" <<
        //"(" << i1+xi <<","<< j1+yi <<","<< k1-1 <<") " <<
        //"(" << i1+xi <<","<< j1+yi <<","<< k1   <<") " <<
        //"(" << i2+xi <<","<< j2+yi <<","<< k2-1 <<") " <<
        //"(" << i2+xi <<","<< j2+yi <<","<< k2   <<") " << "\n";

        atomic_add( yee.jz(i1+xi, j1+yi, k1-1), qvz1*(0.5f - Wzz1[0])*Wxx1[xi+2]*Wyy1[yi+2] );
        atomic_add( yee.jz(i1+xi, j1+yi, k1  ), qvz1*(0.5f + Wzz1[0])*Wxx1[xi+2]*Wyy1[yi+2] );

        atomic_add( yee.jz(i2+xi, j2+yi, k2-1), qvz2*(0.5f - Wzz2[0])*Wxx2[xi+2]*Wyy2[yi+2] );
        atomic_add( yee.jz(i2+xi, j2+yi, k2  ), qvz2*(0.5f + Wzz2[0])*Wxx2[xi+2]*Wyy2[yi+2] );
      }

    }
    //}, con.size(), yee, con);

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
//template class pic::ZigZag<1,3>; // 1D3V
//template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag_2nd<3,3>; // 3D3V

