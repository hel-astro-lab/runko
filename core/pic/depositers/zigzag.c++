#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/pic/depositers/zigzag.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;

template<size_t D, size_t V>
void pic::ZigZag<D,V>::solve( pic::Tile<D>& tile )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& gs = tile.get_grids();
  const auto mins = tile.mins;
  const auto maxs = tile.maxs;

  //clear arrays before new update
  gs.jx.clear();
  gs.jy.clear();
  gs.jz.clear();

  for(auto&& con: tile.containers) {

    const float c = tile.cfl;    // speed of light
    const float q = con.q; // charge

    // skip particle species if zero charge
    if (q == 0.0) continue;
    
    const size_t iy = D >= 2 ? gs.jx.indx(0,1,0) - gs.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? gs.jx.indx(0,0,1) - gs.ex.indx(0,0,0) : 0;

    // no vectorization here since we dont use the general iterator
    //for(size_t n=0; n<con.size(); n++) {

    UniIter::iterate([=] DEVCALLABLE (
                size_t n, 
                emf::Grids &gs,
                pic::ParticleContainer<D>& con
                ){

      //--------------------------------------------------
      float u = con.vel(0,n);
      float v = con.vel(1,n);
      float w = con.vel(2,n);
      float invgam = 1.0/sqrt(1.0 + u*u + v*v + w*w);

      //--------------------------------------------------
      // new (normalized) location, x_{n+1}
      float x2 = D >= 1 ? con.loc(0,n) - mins[0] : 0.0f;
      float y2 = D >= 2 ? con.loc(1,n) - mins[1] : 0.0f;
      float z2 = D >= 3 ? con.loc(2,n) - mins[2] : 0.0f;

      // previos location, x_n
      float x1 = x2 - u*invgam*c;
      float y1 = y2 - v*invgam*c;
      float z1 = z2 - w*invgam*c; 

      //--------------------------------------------------
      int i1  = D >= 1 ? floor(x1) : 0;
      int i2  = D >= 1 ? floor(x2) : 0;
      int j1  = D >= 2 ? floor(y1) : 0;
      int j2  = D >= 2 ? floor(y2) : 0;
      int k1  = D >= 3 ? floor(z1) : 0;
      int k2  = D >= 3 ? floor(z2) : 0;

      // relay point; +1 is equal to +\Delta x
      float xr = min( float(min(i1,i2)+1), max( float(max(i1,i2)), float(0.5f*(x1+x2)) ) );
      float yr = min( float(min(j1,j2)+1), max( float(max(j1,j2)), float(0.5f*(y1+y2)) ) );
      float zr = min( float(min(k1,k2)+1), max( float(max(k1,k2)), float(0.5f*(z1+z2)) ) );

      ////--------------------------------------------------
      //// +q since - sign is already included in the Ampere's equation
      ////q = weight*qe;
      float Fx1 = +q*con.wgt(n)*(xr - x1);
      float Fy1 = +q*con.wgt(n)*(yr - y1);
      float Fz1 = +q*con.wgt(n)*(zr - z1);

      float Fx2 = +q*con.wgt(n)*(x2 - xr);
      float Fy2 = +q*con.wgt(n)*(y2 - yr);
      float Fz2 = +q*con.wgt(n)*(z2 - zr);

      float Wx1 = D >= 1 ? 0.5f*(x1 + xr) - i1 : 0.0f;
      float Wy1 = D >= 2 ? 0.5f*(y1 + yr) - j1 : 0.0f;
      float Wz1 = D >= 3 ? 0.5f*(z1 + zr) - k1 : 0.0f;

      float Wx2 = D >= 1 ? 0.5f*(x2 + xr) - i2 : 0.0f;
      float Wy2 = D >= 2 ? 0.5f*(y2 + yr) - j2 : 0.0f;
      float Wz2 = D >= 3 ? 0.5f*(z2 + zr) - k2 : 0.0f;

     //-------------------------------------------------- 
     // check outflow
#if DEBUG
       
      const int H = 3;
      if((i1 < -H || i1 >= maxs[0] + H-1 ||
          i2 < -H || i2 >= maxs[0] + H-1 ||
          j1 < -H || j1 >= maxs[1] + H-1 ||
          j2 < -H || j2 >= maxs[1] + H-1 ||
          k1 < -H || k1 >= maxs[2] + H-1 ||
          k2 < -H || k2 >= maxs[2] + H-1 ) ){
          //&& con.wgt(n) > 0.0) {

        std::cerr << "ERROR ZIGZAG:" 
                  << " i1 " << i1 << " i2 " << i2
                  << " j1 " << j1 << " j2 " << j2
                  << " k1 " << k1 << " k2 " << k2
                  << " x1 " << x1 << " x2 " << x2
                  << " y1 " << y1 << " y2 " << y2
                  << " z1 " << z1 << " z2 " << z2
                  << " v " << u << " " << v << " " << w 
                  << " wgt " << con.wgt(n) << " info " << con.info(n) 
                  << std::endl;

        // do not deposit anything
        Fx1 = 0.0, Fx2 = 0.0, Fy1 = 0.0, Fy2 = 0.0, Fz1 = 0.0, Fz2 = 0.0;
        i1=0, i2=0, j1=0, j2=0, k1=0, k2=0;
        //assert(false);
      }
#endif

      //--------------------------------------------------
      // one-dimensional indices
        
      const size_t ind1 = gs.jx.indx(i1,j1,k1);
      const size_t ind2 = gs.jx.indx(i2,j2,k2);
        
      if(D>=1) atomic_add( gs.jx(ind1            ), Fx1*(1.0-Wy1)*(1.0-Wz1) );
      if(D>=2) atomic_add( gs.jx(ind1    +iy     ), Fx1*Wy1      *(1.0-Wz1) );
      if(D>=3) atomic_add( gs.jx(ind1        +iz ), Fx1*(1.0-Wy1)*Wz1       );
      if(D>=3) atomic_add( gs.jx(ind1    +iy +iz ), Fx1*Wy1      *Wz1       );

      if(D>=1) atomic_add( gs.jx(ind2            ), Fx2*(1.0-Wy2)*(1.0-Wz2) );
      if(D>=2) atomic_add( gs.jx(ind2    +iy     ), Fx2*Wy2      *(1.0-Wz2) );
      if(D>=3) atomic_add( gs.jx(ind2        +iz ), Fx2*(1.0-Wy2)*Wz2       );
      if(D>=3) atomic_add( gs.jx(ind2    +iy +iz ), Fx2*Wy2      *Wz2       );

      // jy
      if(D>=1) atomic_add( gs.jy(ind1            ), Fy1*(1.0-Wx1)*(1.0-Wz1) );
      if(D>=1) atomic_add( gs.jy(ind1 +1         ), Fy1*Wx1      *(1.0-Wz1) );
      if(D>=3) atomic_add( gs.jy(ind1        +iz ), Fy1*(1.0-Wx1)*Wz1       );
      if(D>=3) atomic_add( gs.jy(ind1 +1     +iz ), Fy1*Wx1      *Wz1       );

      if(D>=1) atomic_add( gs.jy(ind2            ), Fy2*(1.0-Wx2)*(1.0-Wz2) );
      if(D>=1) atomic_add( gs.jy(ind2 +1         ), Fy2*Wx2      *(1.0-Wz2) );
      if(D>=3) atomic_add( gs.jy(ind2        +iz ), Fy2*(1.0-Wx2)*Wz2       );
      if(D>=3) atomic_add( gs.jy(ind2 +1     +iz ), Fy2*Wx2      *Wz2       );

      // jz
      if(D>=1) atomic_add( gs.jz(ind1            ), Fz1*(1.0-Wx1)*(1.0-Wy1) );
      if(D>=1) atomic_add( gs.jz(ind1 +1         ), Fz1*Wx1      *(1.0-Wy1) );
      if(D>=2) atomic_add( gs.jz(ind1    +iy     ), Fz1*(1.0-Wx1)*Wy1       );
      if(D>=2) atomic_add( gs.jz(ind1 +1 +iy     ), Fz1*Wx1      *Wy1       );

      if(D>=1) atomic_add( gs.jz(ind2            ), Fz2*(1.0-Wx2)*(1.0-Wy2) );
      if(D>=1) atomic_add( gs.jz(ind2 +1         ), Fz2*Wx2      *(1.0-Wy2) );
      if(D>=2) atomic_add( gs.jz(ind2    +iy     ), Fz2*(1.0-Wx2)*Wy2       );
      if(D>=2) atomic_add( gs.jz(ind2 +1 +iy     ), Fz2*Wx2      *Wy2       );


      // multid indexing version
      //--------------------------------------------------
      // jx
      //if(D>=1) atomic_add( gs.jx(i1  , j1  , k1  ), Fx1*(1.0-Wy1)*(1.0-Wz1) );
      //if(D>=2) atomic_add( gs.jx(i1  , j1+1, k1  ), Fx1*Wy1      *(1.0-Wz1) );
      //if(D>=3) atomic_add( gs.jx(i1  , j1  , k1+1), Fx1*(1.0-Wy1)*Wz1       );
      //if(D>=3) atomic_add( gs.jx(i1  , j1+1, k1+1), Fx1*Wy1      *Wz1       );

      //if(D>=1) atomic_add( gs.jx(i2  , j2  , k2  ), Fx2*(1.0-Wy2)*(1.0-Wz2) );
      //if(D>=2) atomic_add( gs.jx(i2  , j2+1, k2  ), Fx2*Wy2      *(1.0-Wz2) );
      //if(D>=3) atomic_add( gs.jx(i2  , j2  , k2+1), Fx2*(1.0-Wy2)*Wz2       );
      //if(D>=3) atomic_add( gs.jx(i2  , j2+1, k2+1), Fx2*Wy2      *Wz2       );

      //// jy
      //if(D>=1) atomic_add( gs.jy(i1  , j1  , k1  ), Fy1*(1.0-Wx1)*(1.0-Wz1) );
      //if(D>=1) atomic_add( gs.jy(i1+1, j1  , k1  ), Fy1*Wx1      *(1.0-Wz1) );
      //if(D>=3) atomic_add( gs.jy(i1  , j1  , k1+1), Fy1*(1.0-Wx1)*Wz1       );
      //if(D>=3) atomic_add( gs.jy(i1+1, j1  , k1+1), Fy1*Wx1      *Wz1       );

      //if(D>=1) atomic_add( gs.jy(i2  , j2  , k2  ), Fy2*(1.0-Wx2)*(1.0-Wz2) );
      //if(D>=1) atomic_add( gs.jy(i2+1, j2  , k2  ), Fy2*Wx2      *(1.0-Wz2) );
      //if(D>=3) atomic_add( gs.jy(i2  , j2  , k2+1), Fy2*(1.0-Wx2)*Wz2       );
      //if(D>=3) atomic_add( gs.jy(i2+1, j2  , k2+1), Fy2*Wx2      *Wz2       );

      //// jz
      //if(D>=1) atomic_add( gs.jz(i1  , j1  , k1  ), Fz1*(1.0-Wx1)*(1.0-Wy1) );
      //if(D>=1) atomic_add( gs.jz(i1+1, j1  , k1  ), Fz1*Wx1      *(1.0-Wy1) );
      //if(D>=2) atomic_add( gs.jz(i1  , j1+1, k1  ), Fz1*(1.0-Wx1)*Wy1       );
      //if(D>=2) atomic_add( gs.jz(i1+1, j1+1, k1  ), Fz1*Wx1      *Wy1       );

      //if(D>=1) atomic_add( gs.jz(i2  , j2  , k2  ), Fz2*(1.0-Wx2)*(1.0-Wy2) );
      //if(D>=1) atomic_add( gs.jz(i2+1, j2  , k2  ), Fz2*Wx2      *(1.0-Wy2) );
      //if(D>=2) atomic_add( gs.jz(i2  , j2+1, k2  ), Fz2*(1.0-Wx2)*Wy2       );
      //if(D>=2) atomic_add( gs.jz(i2+1, j2+1, k2  ), Fz2*Wx2      *Wy2       );
      
    }, con.size(), gs, con);


    UniIter::sync();
  }//end of loop over species


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation
template class pic::ZigZag<1,3>; // 1D3V
template class pic::ZigZag<2,3>; // 2D3V
template class pic::ZigZag<3,3>; // 3D3V
