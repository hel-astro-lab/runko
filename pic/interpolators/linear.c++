#include "linear.h"

#include <cmath> 
#include <cassert>
#include <nvtx3/nvToolsExt.h> 
#include "../../tools/iter/iter.h"

template<size_t D, size_t V>
void pic::LinearInterpolator<D,V>::solve(
    pic::Tile<D>& tile)
{
nvtxRangePush(__PRETTY_FUNCTION__);

  // get reference to the Yee grid 
  auto& yee = tile.get_yee();

  for(auto&& container : tile.containers) {

    int nparts = container.size();

    // initialize pointers to particle arrays
    real_prtcl* loc[3];
    for( int i=0; i<3; i++) 
      loc[i] = &( container.loc(i,0) );

    /// resize internal arrays
    container.Epart.resize(3*nparts);
    container.Bpart.resize(3*nparts);
      
    real_prtcl *ex, *ey, *ez, *bx, *by, *bz;
    ex = &( container.Epart[0*nparts] );
    ey = &( container.Epart[1*nparts] );
    ez = &( container.Epart[2*nparts] );

    bx = &( container.Bpart[0*nparts] );
    by = &( container.Bpart[1*nparts] );
    bz = &( container.Bpart[2*nparts] );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    auto mins = tile.mins;
    //auto maxs = tile.maxs;


UniIter::UniIterHost::iterate([=] DEVCALLABLE (int n, fields::YeeLattice &yee){
  //
      int i=0, j=0, k=0;
    real_long dx=0.0, dy=0.0, dz=0.0;
    real_long f,g;

    real_long loc0n, loc1n, loc2n;

    int iz = 1;
    if (D<=2) iz = 0; // flip switch for making array queries 2D

    
      loc0n = static_cast<real_long>( loc[0][n] );
      loc1n = static_cast<real_long>( loc[1][n] );
      loc2n = static_cast<real_long>( loc[2][n] );

      // particle location in the grid
      if (D >= 1) i  = static_cast<int>(floor( loc[0][n] - mins[0] ));
      if (D >= 2) j  = static_cast<int>(floor( loc[1][n] - mins[1] ));
      if (D >= 3) k  = static_cast<int>(floor( loc[2][n] - mins[2] ));

      if (D >= 1) dx = (loc[0][n]-mins[0]) - i;
      if (D >= 2) dy = (loc[1][n]-mins[1]) - j;
      if (D >= 3) dz = (loc[2][n]-mins[2]) - k;

      // TODO: these can be optimized further when we know D
      // TODO: can also be optimized by using 1D indexing (since we can pre-compute index with
      // l = i + iy*(j-1)+iz*(k-1)
      f = yee.ex(i,j,k) + yee.ex(i-1,j,k) +    dx*(yee.ex(i+1,j  ,k   ) - yee.ex(i-1,j  ,k  ));
      f+=                                      dy*(yee.ex(i  ,j+1,k   ) + yee.ex(i-1,j+1,k  )+
                                               dx*(yee.ex(i+1,j+1,k   ) - yee.ex(i-1,j+1,k  ))-f);
      g = yee.ex(i,j,k+1)+yee.ex(i-1,j,k+iz)+  dx*(yee.ex(i+1,j  ,k+iz) - yee.ex(i-1,j  ,k+iz));
      g+=                                      dy*(yee.ex(i  ,j+1,k+iz) + yee.ex(i-1,j+1,k+iz)
                                           +   dx*(yee.ex(i+1,j+1,k+iz) - yee.ex(i-1,j+1,k+iz))-g);
      ex[n] = static_cast<real_prtcl>( 0.5*(f+dz*(g-f)) );

      f = yee.ey(i,j,k)+yee.ey(i,j-1,k)+       dy*(yee.ey(i  ,j+1,k   ) - yee.ey(i  ,j-1,k  ));
      f+=                                      dz*(yee.ey(i  ,j  ,k+iz) + yee.ey(i  ,j-1,k+iz)+      
                                               dy*(yee.ey(i  ,j+1,k+iz) - yee.ey(i  ,j-1,k+iz))-f);
      g = yee.ey(i+1,j,k)+yee.ey(i+1,j-1,k)+   dy*(yee.ey(i+1,j+1,k   ) - yee.ey(i+1,j-1,k   ));
      g+=                                      dz*(yee.ey(i+1,j  ,k+iz) + yee.ey(i+1,j-1,k+iz)+ 
                                               dy*(yee.ey(i+1,j+1,k+iz) - yee.ey(i+1,j-1,k+iz))-g);
      ey[n] = static_cast<real_prtcl>( 0.5*(f+dx*(g-f)) );

      f = yee.ez(i,j,k)+yee.ez(i,j,k-iz)+      dz*(yee.ez(i  ,j  ,k+iz) - yee.ez(i  ,j  ,k-iz));
      f+=                                      dx*(yee.ez(i+1,j  ,k   ) + yee.ez(i+1,j  ,k-iz)+
                                               dz*(yee.ez(i+1,j  ,k+iz) - yee.ez(i+1,j  ,k-iz))-f);
      g = yee.ez(i,j+1,k)+ yee.ez(i,j+1,k-iz)+ dz*(yee.ez(i  ,j+1,k+iz) - yee.ez(i  ,j+1,k-iz));
      g+=                                      dx*(yee.ez(i+1,j+1,k   ) + yee.ez(i+1,j+1,k-iz)+
                                               dz*(yee.ez(i+1,j+1,k+iz) - yee.ez(i+1,j+1,k-iz))-g);
      ez[n] = static_cast<real_prtcl>( 0.5*(f+dy*(g-f)) );

      f = yee.bx(i,j-1,k)  +yee.bx(i,j-1,k-iz )   +dz*(yee.bx(i,j-1,k+iz)   - yee.bx(i,j-1,k-iz));
      f = yee.bx(i,j,k)    +yee.bx(i,j,k-iz)      +dz*(yee.bx(i,j,k+iz)     - yee.bx(i,j,k-iz))+f+dy 
        *(yee.bx(i,j+1,k)  +yee.bx(i,j+1,k-iz)    +dz*(yee.bx(i,j+1,k+iz)   - yee.bx(i,j+1,k-iz))-f);
      g = yee.bx(i+1,j-1,k)+yee.bx(i+1,j-1,k-iz)  +dz*(yee.bx(i+1,j-1,k+iz) - yee.bx(i+1,j-1,k-iz));
      g = yee.bx(i+1,j,k)  +yee.bx(i+1,j,k-iz)    +dz*(yee.bx(i+1,j,k+iz)   - yee.bx(i+1,j,k-iz))
                                           +g     +dy*(yee.bx(i+1,j+1,k)    + yee.bx(i+1,j+1,k-iz)
                                                  +dz*(yee.bx(i+1,j+1,k+iz) - yee.bx(i+1,j+1,k-iz))-g);
      bx[n] = static_cast<real_prtcl>( 0.25*(f+dx*(g-f)) );

      f = yee.by(i,j,k-iz)+yee.by(i-1,j,k-iz)     +dx*(yee.by(i+1,j,k-iz)   - yee.by(i-1,j,k-iz));
      f = yee.by(i,j,k)+yee.by(i-1,j,k)           +dx*(yee.by(i+1,j,k)      - yee.by(i-1,j,k))+f+dz 
        *(yee.by(i,j,k+iz)+yee.by(i-1,j,k+iz)     +dx*(yee.by(i+1,j,k+iz)   - yee.by(i-1,j,k+iz))-f);
      g = yee.by(i,j+1,k-iz)+yee.by(i-1,j+1,k-iz) +dx*(yee.by(i+1,j+1,k-iz) - yee.by(i-1,j+1,k-iz));
      g = yee.by(i,j+1,k)+yee.by(i-1,j+1,k)       +dx*(yee.by(i+1,j+1,k)    - yee.by(i-1,j+1,k))
                                               +g +dz*(yee.by(i,j+1,k+iz)   + yee.by(i-1,j+1,k+iz)
                                                  +dx*(yee.by(i+1,j+1,k+iz) - yee.by(i-1,j+1,k+iz))-g);
      by[n] = static_cast<real_prtcl>( 0.25*(f+dy*(g-f)) );

      f = yee.bz(i-1,j,k)+yee.bz(i-1,j-1,k )      +dy*(yee.bz(i-1,j+1,k)    - yee.bz(i-1,j-1,k));
      f = yee.bz(i,j,k)+yee.bz(i,j-1,k)           +dy*(yee.bz(i,j+1,k)      - yee.bz(i,j-1,k))+f+dx 
        * (yee.bz(i+1,j,k)+yee.bz(i+1,j-1,k)      +dy*(yee.bz(i+1,j+1,k)    - yee.bz(i+1,j-1,k))-f);
      g = yee.bz(i-1,j, k+iz)+yee.bz(i-1,j-1,k+iz)+dy*(yee.bz(i-1,j+1,k+iz) - yee.bz(i-1,j-1,k+iz));
      g = yee.bz(i,j,k+iz)+yee.bz(i,j-1,k+iz )    +dy*(yee.bz(i,j+1,k+iz)   - yee.bz(i,j-1,k+iz))
                                               +g +dx*(yee.bz(i+1,j,k+iz)   + yee.bz(i+1,j-1,k+iz)
                                                  +dy*(yee.bz(i+1,j+1,k+iz) - yee.bz(i+1,j-1,k+iz))-g);
      bz[n] = static_cast<real_prtcl>( 0.25*(f+dz*(g-f)) );
    
}, nparts, yee);
UniIter::sync();
    /*
    // TODO: think SIMD (not possible due to ijk writing to yee) // how is it not possible all writes happen to "n" ?
    #pragma omp parallel for
    for(int n=n1; n<n2; n++) {

    int i=0, j=0, k=0;
    real_long dx=0.0, dy=0.0, dz=0.0;
    real_long f,g;

    real_long loc0n, loc1n, loc2n;

    int iz = 1;
    if (D<=2) iz = 0; // flip switch for making array queries 2D

    
      loc0n = static_cast<real_long>( loc[0][n] );
      loc1n = static_cast<real_long>( loc[1][n] );
      loc2n = static_cast<real_long>( loc[2][n] );

      // particle location in the grid
      if (D >= 1) i  = static_cast<int>(floor( loc[0][n] - mins[0] ));
      if (D >= 2) j  = static_cast<int>(floor( loc[1][n] - mins[1] ));
      if (D >= 3) k  = static_cast<int>(floor( loc[2][n] - mins[2] ));

      if (D >= 1) dx = (loc[0][n]-mins[0]) - i;
      if (D >= 2) dy = (loc[1][n]-mins[1]) - j;
      if (D >= 3) dz = (loc[2][n]-mins[2]) - k;

      // TODO: these can be optimized further when we know D
      // TODO: can also be optimized by using 1D indexing (since we can pre-compute index with
      // l = i + iy*(j-1)+iz*(k-1)
      f = yee.ex(i,j,k) + yee.ex(i-1,j,k) +    dx*(yee.ex(i+1,j  ,k   ) - yee.ex(i-1,j  ,k  ));
      f+=                                      dy*(yee.ex(i  ,j+1,k   ) + yee.ex(i-1,j+1,k  )+
                                               dx*(yee.ex(i+1,j+1,k   ) - yee.ex(i-1,j+1,k  ))-f);
      g = yee.ex(i,j,k+1)+yee.ex(i-1,j,k+iz)+  dx*(yee.ex(i+1,j  ,k+iz) - yee.ex(i-1,j  ,k+iz));
      g+=                                      dy*(yee.ex(i  ,j+1,k+iz) + yee.ex(i-1,j+1,k+iz)
                                           +   dx*(yee.ex(i+1,j+1,k+iz) - yee.ex(i-1,j+1,k+iz))-g);
      ex[n] = static_cast<real_prtcl>( 0.5*(f+dz*(g-f)) );

      f = yee.ey(i,j,k)+yee.ey(i,j-1,k)+       dy*(yee.ey(i  ,j+1,k   ) - yee.ey(i  ,j-1,k  ));
      f+=                                      dz*(yee.ey(i  ,j  ,k+iz) + yee.ey(i  ,j-1,k+iz)+      
                                               dy*(yee.ey(i  ,j+1,k+iz) - yee.ey(i  ,j-1,k+iz))-f);
      g = yee.ey(i+1,j,k)+yee.ey(i+1,j-1,k)+   dy*(yee.ey(i+1,j+1,k   ) - yee.ey(i+1,j-1,k   ));
      g+=                                      dz*(yee.ey(i+1,j  ,k+iz) + yee.ey(i+1,j-1,k+iz)+ 
                                               dy*(yee.ey(i+1,j+1,k+iz) - yee.ey(i+1,j-1,k+iz))-g);
      ey[n] = static_cast<real_prtcl>( 0.5*(f+dx*(g-f)) );

      f = yee.ez(i,j,k)+yee.ez(i,j,k-iz)+      dz*(yee.ez(i  ,j  ,k+iz) - yee.ez(i  ,j  ,k-iz));
      f+=                                      dx*(yee.ez(i+1,j  ,k   ) + yee.ez(i+1,j  ,k-iz)+
                                               dz*(yee.ez(i+1,j  ,k+iz) - yee.ez(i+1,j  ,k-iz))-f);
      g = yee.ez(i,j+1,k)+ yee.ez(i,j+1,k-iz)+ dz*(yee.ez(i  ,j+1,k+iz) - yee.ez(i  ,j+1,k-iz));
      g+=                                      dx*(yee.ez(i+1,j+1,k   ) + yee.ez(i+1,j+1,k-iz)+
                                               dz*(yee.ez(i+1,j+1,k+iz) - yee.ez(i+1,j+1,k-iz))-g);
      ez[n] = static_cast<real_prtcl>( 0.5*(f+dy*(g-f)) );

      f = yee.bx(i,j-1,k)  +yee.bx(i,j-1,k-iz )   +dz*(yee.bx(i,j-1,k+iz)   - yee.bx(i,j-1,k-iz));
      f = yee.bx(i,j,k)    +yee.bx(i,j,k-iz)      +dz*(yee.bx(i,j,k+iz)     - yee.bx(i,j,k-iz))+f+dy 
        *(yee.bx(i,j+1,k)  +yee.bx(i,j+1,k-iz)    +dz*(yee.bx(i,j+1,k+iz)   - yee.bx(i,j+1,k-iz))-f);
      g = yee.bx(i+1,j-1,k)+yee.bx(i+1,j-1,k-iz)  +dz*(yee.bx(i+1,j-1,k+iz) - yee.bx(i+1,j-1,k-iz));
      g = yee.bx(i+1,j,k)  +yee.bx(i+1,j,k-iz)    +dz*(yee.bx(i+1,j,k+iz)   - yee.bx(i+1,j,k-iz))
                                           +g     +dy*(yee.bx(i+1,j+1,k)    + yee.bx(i+1,j+1,k-iz)
                                                  +dz*(yee.bx(i+1,j+1,k+iz) - yee.bx(i+1,j+1,k-iz))-g);
      bx[n] = static_cast<real_prtcl>( 0.25*(f+dx*(g-f)) );

      f = yee.by(i,j,k-iz)+yee.by(i-1,j,k-iz)     +dx*(yee.by(i+1,j,k-iz)   - yee.by(i-1,j,k-iz));
      f = yee.by(i,j,k)+yee.by(i-1,j,k)           +dx*(yee.by(i+1,j,k)      - yee.by(i-1,j,k))+f+dz 
        *(yee.by(i,j,k+iz)+yee.by(i-1,j,k+iz)     +dx*(yee.by(i+1,j,k+iz)   - yee.by(i-1,j,k+iz))-f);
      g = yee.by(i,j+1,k-iz)+yee.by(i-1,j+1,k-iz) +dx*(yee.by(i+1,j+1,k-iz) - yee.by(i-1,j+1,k-iz));
      g = yee.by(i,j+1,k)+yee.by(i-1,j+1,k)       +dx*(yee.by(i+1,j+1,k)    - yee.by(i-1,j+1,k))
                                               +g +dz*(yee.by(i,j+1,k+iz)   + yee.by(i-1,j+1,k+iz)
                                                  +dx*(yee.by(i+1,j+1,k+iz) - yee.by(i-1,j+1,k+iz))-g);
      by[n] = static_cast<real_prtcl>( 0.25*(f+dy*(g-f)) );

      f = yee.bz(i-1,j,k)+yee.bz(i-1,j-1,k )      +dy*(yee.bz(i-1,j+1,k)    - yee.bz(i-1,j-1,k));
      f = yee.bz(i,j,k)+yee.bz(i,j-1,k)           +dy*(yee.bz(i,j+1,k)      - yee.bz(i,j-1,k))+f+dx 
        * (yee.bz(i+1,j,k)+yee.bz(i+1,j-1,k)      +dy*(yee.bz(i+1,j+1,k)    - yee.bz(i+1,j-1,k))-f);
      g = yee.bz(i-1,j, k+iz)+yee.bz(i-1,j-1,k+iz)+dy*(yee.bz(i-1,j+1,k+iz) - yee.bz(i-1,j-1,k+iz));
      g = yee.bz(i,j,k+iz)+yee.bz(i,j-1,k+iz )    +dy*(yee.bz(i,j+1,k+iz)   - yee.bz(i,j-1,k+iz))
                                               +g +dx*(yee.bz(i+1,j,k+iz)   + yee.bz(i+1,j-1,k+iz)
                                                  +dy*(yee.bz(i+1,j+1,k+iz) - yee.bz(i+1,j-1,k+iz))-g);
      bz[n] = static_cast<real_prtcl>( 0.25*(f+dz*(g-f)) );
    }
*/
  } // end of loop over species
nvtxRangePop();
  }


//--------------------------------------------------
// explicit template instantiation

//template class pic::LinearInterpolator<1,3>; // 1D3V
template class pic::LinearInterpolator<2,3>; // 2D3V
template class pic::LinearInterpolator<3,3>; // 3D3V //TODO; validate

