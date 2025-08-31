#include <cmath> 
#include <cassert>

#include "core/pic/interpolators/linear_1st.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif



DEVCALLABLE inline float _lerp(
      float c000, float c100, float c010, float c110,
      float c001, float c101, float c011, float c111,
      float dx, float dy, float dz) 
{
      float c00 = c000 * (1.0-dx) + c100 * dx;
      float c10 = c010 * (1.0-dx) + c110 * dx;
      float c0  = c00  * (1.0-dy) + c10  * dy;
      float c01 = c001 * (1.0-dx) + c101 * dx;
      float c11 = c011 * (1.0-dx) + c111 * dx;
      float c1  = c01  * (1.0-dy) + c11  * dy;
      float c   = c0   * (1.0-dz) + c1   * dz;
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
  auto& gs = tile.get_grids();

  for(auto&& con : tile.containers) {

    /// resize internal arrays
    con.Epart.resize(3*con.size());
    con.Bpart.resize(3*con.size());

    auto mins = tile.mins;

    // mesh sizes for 1D indexing
    const int iy = D >= 2 ? gs.ex.indx(0,1,0) - gs.ex.indx(0,0,0) : 0;
    const int iz = D >= 3 ? gs.ex.indx(0,0,1) - gs.ex.indx(0,0,0) : 0;


    // loop over particles
    UniIter::iterate([=] DEVCALLABLE( 
                size_t n, 
                emf::Grids& gs,
                pic::ParticleContainer<D>& con){
    
      // normalize to tile units
      float loc0n = D >= 1 ? con.loc(0,n) - mins[0] : 0.0f;
      float loc1n = D >= 2 ? con.loc(1,n) - mins[1] : 0.0f;
      float loc2n = D >= 3 ? con.loc(2,n) - mins[2] : 0.0f;

      // particle location in the grid
      int i = (D >= 1) ? floor(loc0n) : 0;
      int j = (D >= 2) ? floor(loc1n) : 0;
      int k = (D >= 3) ? floor(loc2n) : 0;

      float dx = (D >= 1) ? loc0n - i : 0.0f;
      float dy = (D >= 2) ? loc1n - j : 0.0f;
      float dz = (D >= 3) ? loc2n - k : 0.0f;

      // one-dimensional index
      const int ind = gs.ex.indx(i,j,k);
      float c000, c100, c010, c110, c001, c101, c011, c111;

      //ex
      c000 = 0.5*(gs.ex(ind       ) +gs.ex(ind-1      ));
      c100 = 0.5*(gs.ex(ind       ) +gs.ex(ind+1      ));
      c010 = 0.5*(gs.ex(ind+iy    ) +gs.ex(ind-1+iy   ));
      c110 = 0.5*(gs.ex(ind+iy    ) +gs.ex(ind+1+iy   ));
      c001 = 0.5*(gs.ex(ind+iz    ) +gs.ex(ind-1+iz   ));
      c101 = 0.5*(gs.ex(ind+iz    ) +gs.ex(ind+1+iz   ));
      c011 = 0.5*(gs.ex(ind+iy+iz ) +gs.ex(ind-1+iy+iz));
      c111 = 0.5*(gs.ex(ind+iy+iz ) +gs.ex(ind+1+iy+iz));
      con.ex(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      //ey
      c000 = 0.5*(gs.ey(ind       ) +gs.ey(ind-iy     ));
      c100 = 0.5*(gs.ey(ind+1     ) +gs.ey(ind+1-iy   ));
      c010 = 0.5*(gs.ey(ind       ) +gs.ey(ind+iy     ));
      c110 = 0.5*(gs.ey(ind+1     ) +gs.ey(ind+1+iy   ));
      c001 = 0.5*(gs.ey(ind+iz    ) +gs.ey(ind-iy+iz  ));
      c101 = 0.5*(gs.ey(ind+1+iz  ) +gs.ey(ind+1-iy+iz));
      c011 = 0.5*(gs.ey(ind+iz    ) +gs.ey(ind+iy+iz  ));
      c111 = 0.5*(gs.ey(ind+1+iz  ) +gs.ey(ind+1+iy+iz));
      con.ey(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      //ez
      c000 = 0.5*(gs.ez(ind       ) + gs.ez(ind-iz     ));
      c100 = 0.5*(gs.ez(ind+1     ) + gs.ez(ind+1-iz   ));
      c010 = 0.5*(gs.ez(ind+iy    ) + gs.ez(ind+iy-iz  ));
      c110 = 0.5*(gs.ez(ind+1+iy  ) + gs.ez(ind+1+iy-iz));
      c001 = 0.5*(gs.ez(ind       ) + gs.ez(ind+iz     ));
      c101 = 0.5*(gs.ez(ind+1     ) + gs.ez(ind+1+iz   ));
      c011 = 0.5*(gs.ez(ind+iy    ) + gs.ez(ind+iy+iz  ));
      c111 = 0.5*(gs.ez(ind+1+iy  ) + gs.ez(ind+1+iy+iz));
      con.ez(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);


      //-------------------------------------------------- 
      // bx
      c000 = 0.25*( gs.bx(ind)+   gs.bx(ind-iy)+   gs.bx(ind-iz)+      gs.bx(ind-iy-iz));
      c100 = 0.25*( gs.bx(ind+1)+ gs.bx(ind+1-iy)+ gs.bx(ind+1-iz)+    gs.bx(ind+1-iy-iz));
      c001 = 0.25*( gs.bx(ind)+   gs.bx(ind+iz)+   gs.bx(ind-iy)+      gs.bx(ind-iy+iz));
      c101 = 0.25*( gs.bx(ind+1)+ gs.bx(ind+1+iz)+ gs.bx(ind+1-iy)+    gs.bx(ind+1-iy+iz));
      c010 = 0.25*( gs.bx(ind)+   gs.bx(ind+iy)+   gs.bx(ind-iz)+      gs.bx(ind+iy-iz));
      c110 = 0.25*( gs.bx(ind+1)+ gs.bx(ind+1-iz)+ gs.bx(ind+1+iy-iz)+ gs.bx(ind+1+iy));
      c011 = 0.25*( gs.bx(ind)+   gs.bx(ind+iy)+   gs.bx(ind+iy+iz)+   gs.bx(ind+iz));
      c111 = 0.25*( gs.bx(ind+1)+ gs.bx(ind+1+iy)+ gs.bx(ind+1+iy+iz)+ gs.bx(ind+1+iz));
      con.bx(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      // by
      c000 = 0.25*( gs.by(ind-1-iz)+    gs.by(ind-1)+       gs.by(ind-iz)+      gs.by(ind));
      c100 = 0.25*( gs.by(ind-iz)+      gs.by(ind)+         gs.by(ind+1-iz)+    gs.by(ind+1));
      c001 = 0.25*( gs.by(ind-1)+       gs.by(ind-1+iz)+    gs.by(ind)+         gs.by(ind+iz));
      c101 = 0.25*( gs.by(ind)+         gs.by(ind+iz)+      gs.by(ind+1)+       gs.by(ind+1+iz));
      c010 = 0.25*( gs.by(ind-1+iy-iz)+ gs.by(ind-1+iy)+    gs.by(ind+iy-iz)+   gs.by(ind+iy));
      c110 = 0.25*( gs.by(ind+iy-iz)+   gs.by(ind+iy)+      gs.by(ind+1+iy-iz)+ gs.by(ind+1+iy));
      c011 = 0.25*( gs.by(ind-1+iy)+    gs.by(ind-1+iy+iz)+ gs.by(ind+iy)+      gs.by(ind+iy+iz));
      c111 = 0.25*( gs.by(ind+iy)+      gs.by(ind+iy+iz)+   gs.by(ind+1+iy)+    gs.by(ind+1+iy+iz));
      con.by(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

      // bz
      c000 = 0.25*( gs.bz(ind-1-iy)+    gs.bz(ind-1)+       gs.bz(ind-iy)+      gs.bz(ind));
      c100 = 0.25*( gs.bz(ind-iy)+      gs.bz(ind)+         gs.bz(ind+1-iy)+    gs.bz(ind+1));
      c001 = 0.25*( gs.bz(ind-1-iy+iz)+ gs.bz(ind-1+iz)+    gs.bz(ind-iy+iz)+   gs.bz(ind+iz));
      c101 = 0.25*( gs.bz(ind-iy+iz)+   gs.bz(ind+iz)+      gs.bz(ind+1-iy+iz)+ gs.bz(ind+1+iz));
      c010 = 0.25*( gs.bz(ind-1)+       gs.bz(ind-1+iy)+    gs.bz(ind)+         gs.bz(ind+iy));
      c110 = 0.25*( gs.bz(ind)+         gs.bz(ind+iy)+      gs.bz(ind+1)+       gs.bz(ind+1+iy));
      c011 = 0.25*( gs.bz(ind-1+iz)+    gs.bz(ind-1+iy+iz)+ gs.bz(ind+iz)+      gs.bz(ind+iy+iz));
      c111 = 0.25*( gs.bz(ind+iz)+      gs.bz(ind+iy+iz)+   gs.bz(ind+1+iz)+    gs.bz(ind+1+iy+iz));
      con.bz(n) = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

    }, con.size(), gs, con);

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

