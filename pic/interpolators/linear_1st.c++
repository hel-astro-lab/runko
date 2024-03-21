#include "linear_1st.h"

#include <cmath> 
#include <cassert>

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif



DEVCALLABLE inline double _lerp(
      double c000, double c100, double c010, double c110,
      double c001, double c101, double c011, double c111,
      double dx, double dy, double dz) 
{
      double c00 = c000 * (1.0-dx) + c100 * dx;
      double c10 = c010 * (1.0-dx) + c110 * dx;
      double c0  = c00  * (1.0-dy) + c10  * dy;
      double c01 = c001 * (1.0-dx) + c101 * dx;
      double c11 = c011 * (1.0-dx) + c111 * dx;
      double c1  = c01  * (1.0-dy) + c11  * dy;
      double c   = c0   * (1.0-dz) + c1   * dz;
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
  auto& yee = tile.get_grids();

  for(auto&& con : tile.containers) {

    /// resize internal arrays
    con.Epart.resize(3*con.size());
    con.Bpart.resize(3*con.size());

    auto mins = tile.mins;

    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;


    // loop over particles
    UniIter::iterate([=] DEVCALLABLE( 
                size_t n, 
                emf::Grids& yee,
                pic::ParticleContainer<D>& con){

      int i=0, j=0, k=0;
      double dx=0.0, dy=0.0, dz=0.0;
    
      // normalize to tile units
      double loc0n = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      double loc1n = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      double loc2n = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      // particle location in the grid
      if(D >= 1) i = floor(loc0n);
      if(D >= 2) j = floor(loc1n);
      if(D >= 3) k = floor(loc2n);

      if(D >= 1) dx = loc0n - i;
      if(D >= 2) dy = loc1n - j;
      if(D >= 3) dz = loc2n - k;

      // one-dimensional index
      const size_t ind = yee.ex.indx(i,j,k);
      double c000, c100, c010, c110, c001, c101, c011, c111;

      //ex
      c000 = 0.5*(yee.ex(ind       ) +yee.ex(ind-1      ));
      c100 = 0.5*(yee.ex(ind       ) +yee.ex(ind+1      ));
      c010 = 0.5*(yee.ex(ind+iy    ) +yee.ex(ind-1+iy   ));
      c110 = 0.5*(yee.ex(ind+iy    ) +yee.ex(ind+1+iy   ));
      c001 = 0.5*(yee.ex(ind+iz    ) +yee.ex(ind-1+iz   ));
      c101 = 0.5*(yee.ex(ind+iz    ) +yee.ex(ind+1+iz   ));
      c011 = 0.5*(yee.ex(ind+iy+iz ) +yee.ex(ind-1+iy+iz));
      c111 = 0.5*(yee.ex(ind+iy+iz ) +yee.ex(ind+1+iy+iz));
      con.ex(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      //ey
      c000 = 0.5*(yee.ey(ind       ) +yee.ey(ind-iy     ));
      c100 = 0.5*(yee.ey(ind+1     ) +yee.ey(ind+1-iy   ));
      c010 = 0.5*(yee.ey(ind       ) +yee.ey(ind+iy     ));
      c110 = 0.5*(yee.ey(ind+1     ) +yee.ey(ind+1+iy   ));
      c001 = 0.5*(yee.ey(ind+iz    ) +yee.ey(ind-iy+iz  ));
      c101 = 0.5*(yee.ey(ind+1+iz  ) +yee.ey(ind+1-iy+iz));
      c011 = 0.5*(yee.ey(ind+iz    ) +yee.ey(ind+iy+iz  ));
      c111 = 0.5*(yee.ey(ind+1+iz  ) +yee.ey(ind+1+iy+iz));
      con.ey(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      //ez
      c000 = 0.5*(yee.ez(ind       ) + yee.ez(ind-iz     ));
      c100 = 0.5*(yee.ez(ind+1     ) + yee.ez(ind+1-iz   ));
      c010 = 0.5*(yee.ez(ind+iy    ) + yee.ez(ind+iy-iz  ));
      c110 = 0.5*(yee.ez(ind+1+iy  ) + yee.ez(ind+1+iy-iz));
      c001 = 0.5*(yee.ez(ind       ) + yee.ez(ind+iz     ));
      c101 = 0.5*(yee.ez(ind+1     ) + yee.ez(ind+1+iz   ));
      c011 = 0.5*(yee.ez(ind+iy    ) + yee.ez(ind+iy+iz  ));
      c111 = 0.5*(yee.ez(ind+1+iy  ) + yee.ez(ind+1+iy+iz));
      con.ez(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);


      //-------------------------------------------------- 
      // bx
      c000 = 0.25*( yee.bx(ind)+   yee.bx(ind-iy)+   yee.bx(ind-iz)+      yee.bx(ind-iy-iz));
      c100 = 0.25*( yee.bx(ind+1)+ yee.bx(ind+1-iy)+ yee.bx(ind+1-iz)+    yee.bx(ind+1-iy-iz));
      c001 = 0.25*( yee.bx(ind)+   yee.bx(ind+iz)+   yee.bx(ind-iy)+      yee.bx(ind-iy+iz));
      c101 = 0.25*( yee.bx(ind+1)+ yee.bx(ind+1+iz)+ yee.bx(ind+1-iy)+    yee.bx(ind+1-iy+iz));
      c010 = 0.25*( yee.bx(ind)+   yee.bx(ind+iy)+   yee.bx(ind-iz)+      yee.bx(ind+iy-iz));
      c110 = 0.25*( yee.bx(ind+1)+ yee.bx(ind+1-iz)+ yee.bx(ind+1+iy-iz)+ yee.bx(ind+1+iy));
      c011 = 0.25*( yee.bx(ind)+   yee.bx(ind+iy)+   yee.bx(ind+iy+iz)+   yee.bx(ind+iz));
      c111 = 0.25*( yee.bx(ind+1)+ yee.bx(ind+1+iy)+ yee.bx(ind+1+iy+iz)+ yee.bx(ind+1+iz));
      con.bx(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      // by
      c000 = 0.25*( yee.by(ind-1-iz)+    yee.by(ind-1)+       yee.by(ind-iz)+      yee.by(ind));
      c100 = 0.25*( yee.by(ind-iz)+      yee.by(ind)+         yee.by(ind+1-iz)+    yee.by(ind+1));
      c001 = 0.25*( yee.by(ind-1)+       yee.by(ind-1+iz)+    yee.by(ind)+         yee.by(ind+iz));
      c101 = 0.25*( yee.by(ind)+         yee.by(ind+iz)+      yee.by(ind+1)+       yee.by(ind+1+iz));
      c010 = 0.25*( yee.by(ind-1+iy-iz)+ yee.by(ind-1+iy)+    yee.by(ind+iy-iz)+   yee.by(ind+iy));
      c110 = 0.25*( yee.by(ind+iy-iz)+   yee.by(ind+iy)+      yee.by(ind+1+iy-iz)+ yee.by(ind+1+iy));
      c011 = 0.25*( yee.by(ind-1+iy)+    yee.by(ind-1+iy+iz)+ yee.by(ind+iy)+      yee.by(ind+iy+iz));
      c111 = 0.25*( yee.by(ind+iy)+      yee.by(ind+iy+iz)+   yee.by(ind+1+iy)+    yee.by(ind+1+iy+iz));
      con.by(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      // bz
      c000 = 0.25*( yee.bz(ind-1-iy)+    yee.bz(ind-1)+       yee.bz(ind-iy)+      yee.bz(ind));
      c100 = 0.25*( yee.bz(ind-iy)+      yee.bz(ind)+         yee.bz(ind+1-iy)+    yee.bz(ind+1));
      c001 = 0.25*( yee.bz(ind-1-iy+iz)+ yee.bz(ind-1+iz)+    yee.bz(ind-iy+iz)+   yee.bz(ind+iz));
      c101 = 0.25*( yee.bz(ind-iy+iz)+   yee.bz(ind+iz)+      yee.bz(ind+1-iy+iz)+ yee.bz(ind+1+iz));
      c010 = 0.25*( yee.bz(ind-1)+       yee.bz(ind-1+iy)+    yee.bz(ind)+         yee.bz(ind+iy));
      c110 = 0.25*( yee.bz(ind)+         yee.bz(ind+iy)+      yee.bz(ind+1)+       yee.bz(ind+1+iy));
      c011 = 0.25*( yee.bz(ind-1+iz)+    yee.bz(ind-1+iy+iz)+ yee.bz(ind+iz)+      yee.bz(ind+iy+iz));
      c111 = 0.25*( yee.bz(ind+iz)+      yee.bz(ind+iy+iz)+   yee.bz(ind+1+iz)+    yee.bz(ind+1+iy+iz));
      con.bz(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

    }, con.size(), yee, con);

    UniIter::sync();
  } // end of loop over species


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation

template class pic::LinearInterpolator<1,3>; // 1D3V
template class pic::LinearInterpolator<2,3>; // 2D3V
template class pic::LinearInterpolator<3,3>; // 3D3V

