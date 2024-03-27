#include <cmath> 
#include <cassert>

#include "pic/interpolators/cubic_3rd.h"
#include "pic/shapes.h"
#include "tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// 2D specialization
template<>
double pic::CubicInterpolator<2>::compute( 
        double* cx, 
        double* cy, 
        double* /*cz*/, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int j, int /*k*/)
{
  //const size_t ind = f.indx(i,j,0);
  double res = 0.0;

  for( int jl=-1 ; jl<=2 ; jl++ ) {
  for( int il=-1 ; il<=2 ; il++ ) {
    //res += cx[il]*cy[jl] * f(ind + il + (jl*iy) ); // 1D index
    res += cx[il]*cy[jl] * f(i+il, j+jl, 0);
  }}
  return res;
};


// 3D specialization
template<>
double pic::CubicInterpolator<3>::compute( 
        double *cx, 
        double *cy, 
        double *cz, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int j, int k)
{
  //const size_t ind = f.indx(i,j,k);
  double res = 0.0;

  for( int kl=-1 ; kl<=2 ; kl++ ) {
  for( int jl=-1 ; jl<=2 ; jl++ ) {
  for( int il=-1 ; il<=2 ; il++ ) {
    //res += cx[il]*cy[jl]*cz[kl] * f(ind + il + (jl*iy) + (kl*iz) ); // 1D
    res += cx[il]*cy[jl]*cz[kl] * f(i+il, j+jl, k+kl);
  }}}
  return res;
};



template<size_t D>
void pic::CubicInterpolator<D>::solve(
    pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  //NOTE: we need this dummy to enforce template D specialization later on
  //std::array<int, D> dummy; 

  // get reference to the Yee grid 
  auto& gs = tile.get_grids();

  for(auto&& con : tile.containers) {

    /// resize internal arrays
    con.Epart.resize(3*con.size());
    con.Bpart.resize(3*con.size());

    auto mins = tile.mins;

    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? gs.ex.indx(0,1,0) - gs.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? gs.ex.indx(0,0,1) - gs.ex.indx(0,0,0) : 0;


    // loop over particles
    //UniIter::iterate([=] DEVCALLABLE( 
    //            size_t n, 
    //            emf::Grids& gs,
    //            pic::ParticleContainer<D>& con){
    for(size_t n=0; n<con.size(); n++) {
      //--------------------------------------------------
      // indices & locs to grid
      double xpn = (D >= 1) ? con.loc(0, n) - mins[0] : con.loc(0, n) ;
      double ypn = (D >= 2) ? con.loc(1, n) - mins[1] : con.loc(1, n) ;
      double zpn = (D >= 3) ? con.loc(2, n) - mins[2] : con.loc(2, n) ;

      //double u = con.vel(0,n);
      //double v = con.vel(1,n);
      //double w = con.vel(2,n);

      // particle location in the primary and dual 1/2-shifted grid
      int ip = floor(xpn)-1.0;
      int jp = floor(ypn)-1.0;
      int kp = floor(zpn)-1.0;
      int id = floor(xpn);
      int jd = floor(ypn);
      int kd = floor(zpn);

      //--------------------------------------------------
      // coefficients on both prime and dual (staggered +0.5) grids
      double cxd[5] = {0.0}, 
             cxp[5] = {0.0}, 
             cyd[5] = {0.0}, 
             cyp[5] = {0.0}, 
             czd[5] = {0.0}, 
             czp[5] = {0.0};

      
      // \Delta x from primary and staggered grid points
      double dxp = xpn-ip;
      double dyp = ypn-jp;
      double dzp = zpn-kp;

      double dxd = xpn-id-0.5;
      double dyd = ypn-jd-0.5;
      double dzd = zpn-kd-0.5;
        
      //--------------------------------------------------
      // compute shape function weights
        
      //ver2: Eneryg conserving Sokolov alternating shape function scheme
      if(D >= 1) W2nd(dxp, &cxp[0] );
      if(D >= 2) W2nd(dyp, &cyp[0] );
      if(D >= 3) W2nd(dzp, &czp[0] );
      if(D >= 1) W3rd(dxd, &cxd[0] );
      if(D >= 2) W3rd(dyd, &cyd[0] );
      if(D >= 3) W3rd(dzd, &czd[0] );

      // default scheme
      //if(D >= 1) W3rd(dxp, &cxp[0] );
      //if(D >= 2) W3rd(dyp, &cyp[0] );
      //if(D >= 3) W3rd(dzp, &czp[0] );
      //if(D >= 1) W3rd(dxd, &cxd[0] );
      //if(D >= 2) W3rd(dyd, &cyd[0] );
      //if(D >= 3) W3rd(dzd, &czd[0] );


      //--------------------------------------------------
      // integrate over shape function, i.e. interpolate
      // NOTE: we shift coefficient array pointers +1 to take into account -1,0,+1 access pattern
      con.ex(n) = compute( &cxd[2], &cyp[2], &czp[2], gs.ex, iy,iz,  id,jp,kp); // Ex(d,p,p)
      con.ey(n) = compute( &cxp[2], &cyd[2], &czp[2], gs.ey, iy,iz,  ip,jd,kp); // Ey(p,d,p)
      con.ez(n) = compute( &cxp[2], &cyp[2], &czd[2], gs.ez, iy,iz,  ip,jp,kd); // Ez(p,p,d)
      con.bx(n) = compute( &cxp[2], &cyd[2], &czd[2], gs.bx, iy,iz,  ip,jd,kd); // Bx(p,d,d)
      con.by(n) = compute( &cxd[2], &cyp[2], &czd[2], gs.by, iy,iz,  id,jp,kd); // By(d,p,d)
      con.bz(n) = compute( &cxd[2], &cyd[2], &czp[2], gs.bz, iy,iz,  id,jd,kp); // Bz(d,d,p)

    //}, con.size(), gs, con);
    }

    UniIter::sync();
  } // end of loop over species


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation

template class pic::CubicInterpolator<2>; // 2D3V
template class pic::CubicInterpolator<3>; // 3D3V
