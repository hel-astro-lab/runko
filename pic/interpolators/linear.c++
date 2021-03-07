#include "linear.h"

#include <cmath> 
#include <cassert>

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


inline real_long _lerp(
      real_long c000,
      real_long c100,
      real_long c010,
      real_long c110,
      real_long c001,
      real_long c101,
      real_long c011,
      real_long c111,
      real_long dx, real_long dy, real_long dz
      ) 
{
      real_long c00 = c000 * (1.0-dx) + c100 * dx;
      real_long c10 = c010 * (1.0-dx) + c110 * dx;
      real_long c0  = c00  * (1.0-dy) + c10  * dy;
      real_long c01 = c001 * (1.0-dx) + c101 * dx;
      real_long c11 = c011 * (1.0-dx) + c111 * dx;
      real_long c1  = c01  * (1.0-dy) + c11  * dy;
      real_long c   = c0   * (1.0-dz) + c1   * dz;
      return c;
}


template<size_t D, size_t V>
void pic::LinearInterpolator<D,V>::solve(
    pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  // get reference to the Yee grid 
  auto& yee = tile.get_yee();

  auto& ex = yee.ex;
  auto& ey = yee.ey;
  auto& ez = yee.ez;

  auto& bx = yee.bx;
  auto& by = yee.by;
  auto& bz = yee.bz;


  for(auto&& container : tile.containers) {

    int nparts = container.size();

    // initialize pointers to particle arrays
    real_prtcl* loc[3];
    for( int i=0; i<3; i++) 
      loc[i] = &( container.loc(i,0) );

    /// resize internal arrays
    container.Epart.resize(3*nparts);
    container.Bpart.resize(3*nparts);
      
    real_prtcl *exn, *eyn, *ezn, *bxn, *byn, *bzn;
    exn = &( container.Epart[0*nparts] );
    eyn = &( container.Epart[1*nparts] );
    ezn = &( container.Epart[2*nparts] );

    bxn = &( container.Bpart[0*nparts] );
    byn = &( container.Bpart[1*nparts] );
    bzn = &( container.Bpart[2*nparts] );

    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    auto mins = tile.mins;
    //auto maxs = tile.maxs;

    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;


    // loop over particles
    UniIter::iterate([=] DEVCALLABLE (int n, fields::YeeLattice &yee){
      int i=0, j=0, k=0;
      real_long dx=0.0, dy=0.0, dz=0.0;
    
      real_long loc0n, loc1n, loc2n;
      loc0n = static_cast<real_long>( loc[0][n] );
      loc1n = static_cast<real_long>( loc[1][n] );
      loc2n = static_cast<real_long>( loc[2][n] );

      // particle location in the grid
      // NOTE: trunc() ensures that prtcls outside the tile do not crash this loop 
      // (because it rounds e.g. xloc=-0.1 => i=0 and dx=-0.1). They should be
      // automatically cleaned on next time step to their real tiles.
      if(D >= 1) i  = static_cast<int>(floor(loc0n));
      if(D >= 2) j  = static_cast<int>(floor(loc1n));
      if(D >= 3) k  = static_cast<int>(floor(loc2n));

      if(D >= 1) dx = loc0n - i;
      if(D >= 2) dy = loc1n - j;
      if(D >= 3) dz = loc2n - k;

      // normalize to tile units
      if(D >= 1) i -= mins[0];
      if(D >= 2) j -= mins[1];
      if(D >= 3) k -= mins[2];

      // one-dimensional index
      const size_t ind = yee.ex.indx(i,j,k);

      real_long c000, c100, c010, c110, c001, c101, c011, c111;

      //ex
      c000 = 0.5*(ex(ind       ) +ex(ind-1      ));
      c100 = 0.5*(ex(ind       ) +ex(ind+1      ));
      c010 = 0.5*(ex(ind+iy    ) +ex(ind-1+iy   ));
      c110 = 0.5*(ex(ind+iy    ) +ex(ind+1+iy   ));
      c001 = 0.5*(ex(ind+iz    ) +ex(ind-1+iz   ));
      c101 = 0.5*(ex(ind+iz    ) +ex(ind+1+iz   ));
      c011 = 0.5*(ex(ind+iy+iz ) +ex(ind-1+iy+iz));
      c111 = 0.5*(ex(ind+iy+iz ) +ex(ind+1+iy+iz));
      exn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );

      //ey
      c000 = 0.5*(ey(ind      ) +ey(ind-iy     ));
      c100 = 0.5*(ey(ind+1    ) +ey(ind+1-iy   ));
      c010 = 0.5*(ey(ind      ) +ey(ind+iy     ));
      c110 = 0.5*(ey(ind+1    ) +ey(ind+1+iy   ));
      c001 = 0.5*(ey(ind+iz   ) +ey(ind-iy+iz  ));
      c101 = 0.5*(ey(ind+1+iz ) +ey(ind+1-iy+iz));
      c011 = 0.5*(ey(ind+iz   ) +ey(ind+iy+iz  ));
      c111 = 0.5*(ey(ind+1+iz ) +ey(ind+1+iy+iz));
      eyn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );

      //ez
      c000 = 0.5*(ez(ind      ) + ez(ind-iz     ));
      c100 = 0.5*(ez(ind+1    ) + ez(ind+1-iz   ));
      c010 = 0.5*(ez(ind+iy   ) + ez(ind+iy-iz  ));
      c110 = 0.5*(ez(ind+1+iy ) + ez(ind+1+iy-iz));
      c001 = 0.5*(ez(ind      ) + ez(ind+iz     ));
      c101 = 0.5*(ez(ind+1    ) + ez(ind+1+iz   ));
      c011 = 0.5*(ez(ind+iy   ) + ez(ind+iy+iz  ));
      c111 = 0.5*(ez(ind+1+iy ) + ez(ind+1+iy+iz));
      ezn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );


      //-------------------------------------------------- 

      // bx
      c000 = 0.25*( bx(ind)+   bx(ind-iy)+   bx(ind-iz)+      bx(ind-iy-iz));
      c100 = 0.25*( bx(ind+1)+ bx(ind+1-iy)+ bx(ind+1-iz)+    bx(ind+1-iy-iz));
      c001 = 0.25*( bx(ind)+   bx(ind+iz)+   bx(ind-iy)+      bx(ind-iy+iz));
      c101 = 0.25*( bx(ind+1)+ bx(ind+1+iz)+ bx(ind+1-iy)+    bx(ind+1-iy+iz));
      c010 = 0.25*( bx(ind)+   bx(ind+iy)+   bx(ind-iz)+      bx(ind+iy-iz));
      c110 = 0.25*( bx(ind+1)+ bx(ind+1-iz)+ bx(ind+1+iy-iz)+ bx(ind+1+iy));
      c011 = 0.25*( bx(ind)+   bx(ind+iy)+   bx(ind+iy+iz)+   bx(ind+iz));
      c111 = 0.25*( bx(ind+1)+ bx(ind+1+iy)+ bx(ind+1+iy+iz)+ bx(ind+1+iz));
      bxn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );

      // by
      c000 = 0.25*( by(ind-1-iz)+    by(ind-1)+       by(ind-iz)+      by(ind));
      c100 = 0.25*( by(ind-iz)+      by(ind)+         by(ind+1-iz)+    by(ind+1));
      c001 = 0.25*( by(ind-1)+       by(ind-1+iz)+    by(ind)+         by(ind+iz));
      c101 = 0.25*( by(ind)+         by(ind+iz)+      by(ind+1)+       by(ind+1+iz));
      c010 = 0.25*( by(ind-1+iy-iz)+ by(ind-1+iy)+    by(ind+iy-iz)+   by(ind+iy));
      c110 = 0.25*( by(ind+iy-iz)+   by(ind+iy)+      by(ind+1+iy-iz)+ by(ind+1+iy));
      c011 = 0.25*( by(ind-1+iy)+    by(ind-1+iy+iz)+ by(ind+iy)+      by(ind+iy+iz));
      c111 = 0.25*( by(ind+iy)+      by(ind+iy+iz)+   by(ind+1+iy)+    by(ind+1+iy+iz));
      byn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );

      // bz
      c000 = 0.25*( bz(ind-1-iy)+    bz(ind-1)+       bz(ind-iy)+      bz(ind));
      c100 = 0.25*( bz(ind-iy)+      bz(ind)+         bz(ind+1-iy)+    bz(ind+1));
      c001 = 0.25*( bz(ind-1-iy+iz)+ bz(ind-1+iz)+    bz(ind-iy+iz)+   bz(ind+iz));
      c101 = 0.25*( bz(ind-iy+iz)+   bz(ind+iz)+      bz(ind+1-iy+iz)+ bz(ind+1+iz));
      c010 = 0.25*( bz(ind-1)+       bz(ind-1+iy)+    bz(ind)+         bz(ind+iy));
      c110 = 0.25*( bz(ind)+         bz(ind+iy)+      bz(ind+1)+       bz(ind+1+iy));
      c011 = 0.25*( bz(ind-1+iz)+    bz(ind-1+iy+iz)+ bz(ind+iz)+      bz(ind+iy+iz));
      c111 = 0.25*( bz(ind+iz)+      bz(ind+iy+iz)+   bz(ind+1+iz)+    bz(ind+1+iy+iz));
      bzn[n] = static_cast<real_prtcl>( _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz) );
    
    }, nparts, yee);

#ifdef GPU
    UniIter::sync();
#endif

  } // end of loop over species

#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation

//template class pic::LinearInterpolator<1,3>; // 1D3V
template class pic::LinearInterpolator<2,3>; // 2D3V
template class pic::LinearInterpolator<3,3>; // 3D3V //TODO; validate

