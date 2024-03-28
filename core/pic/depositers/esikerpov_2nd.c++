#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/pic/depositers/esikerpov_2nd.h"
#include "core/pic/shapes.h"
#include "external/iter/iter.h"


#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using std::min;
using std::max;


template<size_t D, size_t V>
void pic::Esikerpov_2nd<D,V>::solve( pic::Tile<D>& tile )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  auto& gs = tile.get_grids();
  const auto mins = tile.mins;

  //clear arrays before new update
  gs.jx.clear();
  gs.jy.clear();
  gs.jz.clear();

  for(auto&& con : tile.containers) {

    const double c = tile.cfl;    // speed of light
    const double q = con.q; // charge
                            //
    // skip particle species if zero charge
    if (q == 0.0) continue;

    //UniIter::iterate([=] DEVCALLABLE (
    //            size_t n, 
    //            emf::Grids &gs,
    //            pic::ParticleContainer<D>& con
    //            ){
    for(size_t n=0; n<con.size(); n++) {

      // shape arrays
      double 
          Sx1[5] = {0}, Sx2[5] = {0}, 
          Sy1[5] = {0}, Sy2[5] = {0}, 
          Sz1[5] = {0}, Sz2[5] = {0}, 
          DSx[5] = {0}, 
          DSy[5] = {0}, 
          DSz[5] = {0}; 

      // temporary helper arrays
      double tmpJx[5][5], tmpJy[5][5], tmpJz[5][5];
      for(int i=0; i < 5; i++) {
        for(int j=0; j < 5; j++) {
            tmpJx[i][j] = 0.0;
            tmpJy[i][j] = 0.0;
            tmpJz[i][j] = 0.0;
        }
      }


      //--------------------------------------------------
      // particle locations
        
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
      int i2p = D >= 1 ? round(x2) : 0;

      int j1p = D >= 2 ? round(y1) : 0;
      int j2p = D >= 2 ? round(y2) : 0;

      int k1p = D >= 3 ? round(z1) : 0;
      int k2p = D >= 3 ? round(z2) : 0;

      //--------------------------------------------------
      // Esikerpov weights at old position; 
      // NOTE: we center the array to S[2] = W_i by usign +1 offset
      //
      // Argumetn is \Delta x
      // TODO: add if(D >= X)
      W2nd( x1 - double(i1p), &Sx1[1]); 
      W2nd( y1 - double(j1p), &Sy1[1]); 
      W2nd( z1 - double(k1p), &Sz1[1]); 

      // Esikerpov weights at new position
      // NTOE: ooffset is selected based on \Delta I = change of grid cell
      // and then aligned to S[2] with the +1 offset
      // NOTE: these are only needed for the calculation of DS 
      W2nd( x2 - double(i2p), &Sx2[ i2p - i1p + 1 ]); 
      W2nd( y2 - double(j2p), &Sy2[ j2p - j1p + 1 ]); 
      W2nd( z2 - double(k2p), &Sz2[ k2p - k1p + 1 ]); 

      // compute differences
      for( unsigned int i=0; i < 5; i++) {
        DSx[i] = Sx2[i] - Sx1[i];
        DSy[i] = Sy2[i] - Sy1[i];
        DSz[i] = Sz2[i] - Sz1[i];
      }
    
      //-------------------------------------------------- 

      // current calculation
      int iloc, jloc, kloc; //, linindex;
      const int offset = 2; // -2 comes from 2nd order scheme

      // Jx^(d,p,p)
      for(int k=0 ; k<5 ; k++){
        kloc = k + k1p - offset;
        for(int j=0 ; j<5 ; j++){
          jloc = j + j1p - offset;
          for(int i=1 ; i<5 ; i++){
            iloc = i + i1p - offset; 

            tmpJx[j][k] -= q*c * DSx[i-1]*( Sy1[j] * Sz1[k] 
                           + DSy[j]*Sz1[k]/2.
                           + DSz[k]*Sy1[j]/2.
                           + DSy[j]*DSz[k]/3.);

            // TODO 1d indexing
            atomic_add( gs.jx(iloc, jloc, kloc), tmpJx[j][k] );
            //atomic_add( gs.jx(i2p+i-offset, j2p+j-offset, k2p+k-offset), tmpJx[j][k] );
      }}}

      //-------------------------------------------------- 
      // Jy^(p,d,p)
      for(int k=0 ; k<5 ; k++){
        kloc = k + k1p - offset;
        for(int j=1 ; j<5 ; j++){
          jloc = j + j1p - offset;
          for(int i=0 ; i<5 ; i++){
            iloc = i + i1p - offset; 

            tmpJy[i][k] -= q*c * DSy[j-1] * ( Sz1[k]*Sx1[i] 
                           + DSz[k]*Sx1[i]/2.
                           + DSx[i]*Sz1[k]/2.
                           + DSz[k]*DSx[i]/3.);


            atomic_add( gs.jy(iloc, jloc, kloc), tmpJy[i][k] );
            //atomic_add( gs.jy(i2p+i-offset, j2p+j-offset, k2p+k-offset), tmpJy[i][k] );
      }}}


      //-------------------------------------------------- 
      // Jz^(p,p,d)
      for(int k=1 ; k<5 ; k++){
        kloc = k + k1p - offset;
        for(int j=0 ; j<5 ; j++){
          jloc = j + j1p - offset;
          for(int i=0 ; i<5 ; i++){
            iloc = i + i1p - offset; 

            tmpJz[i][j] -= q*c * DSz[k-1] * ( Sx1[i]*Sy1[j] 
                           + DSx[i]*Sy1[j]/2. 
                           + DSy[j]*Sx1[i]/2. 
                           + DSx[i]*DSy[j]/3.);


            atomic_add( gs.jz(iloc, jloc, kloc), tmpJz[i][j] );
            //atomic_add( gs.jz(i2p+i-offset, j2p+j-offset, k2p+k-offset), tmpJz[i][j] );
      }}}



      //-------------------------------------------------- 
      //-------------------------------------------------- 

      bool debug = false;
      for(int iii=0; iii<5; iii++){
          if(Sx1[iii] < 0.0)  debug = true;
          if(Sy1[iii] < 0.0)  debug = true;
          if(Sz1[iii] < 0.0)  debug = true;
          if(Sx2[iii] < 0.0)  debug = true;
          if(Sy2[iii] < 0.0)  debug = true;
          if(Sz2[iii] < 0.0)  debug = true;
      }

      if(debug){
        std::cout 
          << "cur " 
          << " i1: "  << "(" << i1p << "," << j1p << "," << k1p << ")"
          << " i2: "  << "(" << i2p << "," << j2p << "," << k2p << ")"
          << " x1: " << "(" << x1 << "," << y1 << "," << z1 << ")"
          << " x2: " << "(" << x2 << "," << y2 << "," << z2 << ")"
          //<< " dx1d: "<< "(" << dx1d << "," << dy1d << "," << dz1d << ")"
          //<< " dx2d: "<< "(" << dx2d << "," << dy2d << "," << dz2d << ")"
          << " W1x1: "<< "(" << Sx1[0] << "," << Sx1[1] << "," << Sx1[2] << "," << Sx1[3] << "," << Sx1[4] << ")"
          << " W1y1: "<< "(" << Sy1[0] << "," << Sy1[1] << "," << Sy1[2] << "," << Sy1[3] << "," << Sy1[4] << ")"
          << " W1z1: "<< "(" << Sz1[0] << "," << Sz1[1] << "," << Sz1[2] << "," << Sz1[3] << "," << Sz1[4] << ")"
          << " W2x1: "<< "(" << Sx2[0] << "," << Sx2[1] << "," << Sx2[2] << "," << Sx1[3] << "," << Sx1[4] << ")"
          << " W2y1: "<< "(" << Sy2[0] << "," << Sy2[1] << "," << Sy2[2] << "," << Sy1[3] << "," << Sy1[4] << ")"
          << " W2z1: "<< "(" << Sz2[0] << "," << Sz2[1] << "," << Sz2[2] << "," << Sz1[3] << "," << Sz1[4] << ")"
          << "\n";

      }


    }
    //}, con.size(), gs, con);

  }//end of loop over species

}


//--------------------------------------------------
// explicit template instantiation
//template class pic::Esikerpov_2nd<1,3>; // 1D3V
//template class pic::Esikerpov_2nd<2,3>; // 2D3V
template class pic::Esikerpov_2nd<3,3>; // 3D3V

