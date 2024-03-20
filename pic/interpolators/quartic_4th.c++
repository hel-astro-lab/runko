#include "quartic_4th.h"

#include <cmath> 
#include <cassert>

#include "../shapes.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// 2D specialization
template<>
double pic::QuarticInterpolator<2>::compute( 
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

  for( int jl=-2 ; jl<=2 ; jl++ ) {
  for( int il=-2 ; il<=2 ; il++ ) {
    //res += cx[il]*cy[jl] * f(ind + il + (jl*iy) ); // 1D index
    res += cx[il]*cy[jl] * f(i+il, j+jl, 0);
  }}
  return res;
};


// 3D specialization
template<>
double pic::QuarticInterpolator<3>::compute( 
        double* cx, 
        double* cy, 
        double* cz, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int j, int k)
{
  //const size_t ind = f.indx(i,j,k);
  double res = 0.0;

  for(int kl=-2; kl<=2; kl++ ) {
  for(int jl=-2; jl<=2; jl++ ) {
  for(int il=-2; il<=2; il++ ) {
    //res += cx[il]*cy[jl]*cz[kl] * f(ind + il + (jl*iy) + (kl*iz) ); // 1D
    res += cx[il]*cy[jl]*cz[kl] * f(i+il, j+jl, k+kl);
    //res += *(cx+il)* *(cy+jl) * *(cz+kl) * f(i+il, j+jl, k+kl);
  }}}
  return res;
};



template<size_t D>
void pic::QuarticInterpolator<D>::solve(
    pic::Tile<D>& tile)
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  //NOTE: we need this dummy to enforce template D specialization later on
  //std::array<int, D> dummy; 

  // get reference to the Yee grid 
  auto& yee = tile.get_yee();

  for(auto&& con : tile.containers) {

    /// resize internal arrays
    con.Epart.resize(3*con.size());
    con.Bpart.resize(3*con.size());

    auto mins = tile.mins;

    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;


    // loop over particles
    //UniIter::iterate([=] DEVCALLABLE( 
    //            size_t n, 
    //            emf::YeeLattice& yee,
    //            pic::ParticleContainer<D>& con){
    for(size_t n=0; n<con.size(); n++) {
      //--------------------------------------------------
      // indices & locs to grid
      double xpn = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      double ypn = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      double zpn = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      //double u = con.vel(0,n);
      //double v = con.vel(1,n);
      //double w = con.vel(2,n);

      // particle location in the primary and dual 1/2-shifted grid

      // Default scheme
      int ip = round(xpn);
      int jp = round(ypn);
      int kp = round(zpn);
      int id = round(xpn-0.5);
      int jd = round(ypn-0.5);
      int kd = round(zpn-0.5);

      // Sokolev scheme
      //int ip = round(xpn);
      //int jp = round(ypn);
      //int kp = round(zpn);
      //int id = floor(xpn-0.5);
      //int jd = floor(ypn-0.5);
      //int kd = floor(zpn-0.5);

      //--------------------------------------------------
      // coefficients on both prime and dual (staggered +0.5) grids
      double cxd[5] = {0.0}, 
             cxp[5] = {0.0}, 
             cyd[5] = {0.0}, 
             cyp[5] = {0.0}, 
             czd[5] = {0.0}, 
             czp[5] = {0.0};
      
      // \Delta x from primary and staggered grid points
        
      // Default scheme
      double dxp = xpn-ip;
      double dyp = ypn-jp;
      double dzp = zpn-kp;
      double dxd = xpn-id-0.5;
      double dyd = ypn-jd-0.5;
      double dzd = zpn-kd-0.5;

      //--------------------------------------------------
      // compute shape function weights
        
      // default scheme
      if(D >= 1) W4th(dxp, &cxp[0] );
      if(D >= 2) W4th(dyp, &cyp[0] );
      if(D >= 3) W4th(dzp, &czp[0] );
      if(D >= 1) W4th(dxd, &cxd[0] );
      if(D >= 2) W4th(dyd, &cyd[0] );
      if(D >= 3) W4th(dzd, &czd[0] );
      
      //ver2: Eneryg conserving Sokolov alternating shape function scheme
      //if(D >= 1) W3rd(dxd, &cxd[0] );
      //if(D >= 2) W3rd(dyd, &cyd[0] );
      //if(D >= 3) W3rd(dzd, &czd[0] );
      //if(D >= 1) W4th(dxp, &cxp[0] );
      //if(D >= 2) W4th(dyp, &cyp[0] );
      //if(D >= 3) W4th(dzp, &czp[0] );

      bool debug = false;
      for(int iii=0; iii<5; iii++){
          if(cxp[iii] < -1.0e-6)  debug = true;
          if(cxd[iii] < -1.0e-6)  debug = true;
                          
          if(cyp[iii] < -1.0e-6)  debug = true;
          if(cyd[iii] < -1.0e-6)  debug = true;
                          
          if(czp[iii] < -1.0e-6)  debug = true;
          if(czd[iii] < -1.0e-6)  debug = true;

          if(cxp[iii] > 1.0)  debug = true;
          if(cxd[iii] > 1.0)  debug = true;
          if(cyp[iii] > 1.0)  debug = true;
          if(cyd[iii] > 1.0)  debug = true;
          if(czp[iii] > 1.0)  debug = true;
          if(czd[iii] > 1.0)  debug = true;
      }

      if(debug){
        std::cout 
          << "\n\ninterp xyz: " << "(" << xpn << "," << ypn << "," << zpn << ")"
          << " ip: "  << "(" << ip << "," << jp << "," << kp << ")"
          << " id: "  << "(" << id << "," << jd << "," << kd << ")"
          << " dxp: " << "(" << xpn - ip << "," << ypn - jp << "," << zpn - kp << ")"
          << " dxd: " << "(" << xpn - id + 0.5f << "," << ypn - jd + 0.5f << "," << zpn - kd + 0.5f << ")"
          << " W1x: " << "(" << cxp[0] << "," << cxp[1] << "," << cxp[2] << "," << cxp[3] << "," << cxp[4] << ")"
          << " W2x: " << "(" << cxd[0] << "," << cxd[1] << "," << cxd[2] << "," << cxd[3] << "," << cxd[4] << ")"
          << " W1y: " << "(" << cyp[0] << "," << cyp[1] << "," << cyp[2] << "," << cyp[3] << "," << cyp[4] << ")"
          << " W2y: " << "(" << cyd[0] << "," << cyd[1] << "," << cyd[2] << "," << cyd[3] << "," << cyd[4] << ")"
          << " W1z: " << "(" << czp[0] << "," << czp[1] << "," << czp[2] << "," << czp[3] << "," << czp[4] << ")"
          << " W2z: " << "(" << czd[0] << "," << czd[1] << "," << czd[2] << "," << czd[3] << "," << czd[4] << ")"
          << "\n";
      }

      //--------------------------------------------------
      // integrate over shape function, i.e. interpolate
      // NOTE: we shift coefficient array pointers +1 to take into account 
      // 3rd order -1,0,+1,+2 access pattern
      // 4th order -2,-1,0,1,2 access pattern
        
      // default scheme
      con.ex(n) = compute( &cxd[2], &cyp[2], &czp[2], yee.ex, iy,iz,  id,jp,kp); // Ex(d,p,p)
      con.ey(n) = compute( &cxp[2], &cyd[2], &czp[2], yee.ey, iy,iz,  ip,jd,kp); // Ey(p,d,p)
      con.ez(n) = compute( &cxp[2], &cyp[2], &czd[2], yee.ez, iy,iz,  ip,jp,kd); // Ez(p,p,d)
      con.bx(n) = compute( &cxp[2], &cyd[2], &czd[2], yee.bx, iy,iz,  ip,jd,kd); // Bx(p,d,d)
      con.by(n) = compute( &cxd[2], &cyp[2], &czd[2], yee.by, iy,iz,  id,jp,kd); // By(d,p,d)
      con.bz(n) = compute( &cxd[2], &cyd[2], &czp[2], yee.bz, iy,iz,  id,jd,kp); // Bz(d,d,p)

      // Sokolov
      //con.ex(n) = compute( &cxd[2], &cyp[2], &czp[2], yee.ex, iy,iz,  id,jp,kp); // Ex(d,p,p)
      //con.ey(n) = compute( &cxp[2], &cyd[2], &czp[2], yee.ey, iy,iz,  ip,jd,kp); // Ey(p,d,p)
      //con.ez(n) = compute( &cxp[2], &cyp[2], &czd[2], yee.ez, iy,iz,  ip,jp,kd); // Ez(p,p,d)
      //con.bx(n) = compute( &cxp[2], &cyd[2], &czd[2], yee.bx, iy,iz,  ip,jd,kd); // Bx(p,d,d)
      //con.by(n) = compute( &cxd[2], &cyp[2], &czd[2], yee.by, iy,iz,  id,jp,kd); // By(d,p,d)
      //con.bz(n) = compute( &cxd[2], &cyd[2], &czp[2], yee.bz, iy,iz,  id,jd,kp); // Bz(d,d,p)

    //}, con.size(), yee, con);
    }

    UniIter::sync();
  } // end of loop over species


#ifdef GPU
  nvtxRangePop();
#endif

}


//--------------------------------------------------
// explicit template instantiation

template class pic::QuarticInterpolator<2>; // 2D3V
template class pic::QuarticInterpolator<3>; // 3D3V
