#include <cmath> 
#include <cassert>

#include "core/pic/interpolators/quadratic_2nd.h"
#include "core/pic/shapes.h"
#include "external/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// 1D specialization
template<>
double pic::QuadraticInterpolator<1>::compute( 
        double* cx, 
        double* /*cy*/, 
        double* /*cz*/, 
        const toolbox::Mesh<float, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int /*j*/, int /*k*/)
{
  double res = 0.0;
  for( int il=-1 ; il<=1 ; il++ ) {
    res += cx[il] * f(i+il, 0, 0);
  }
  return res;
};


// 2D specialization
template<>
double pic::QuadraticInterpolator<2>::compute( 
        double* cx, 
        double* cy, 
        double* /*cz*/, 
        const toolbox::Mesh<float, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int j, int /*k*/)
{
  //const size_t ind = f.indx(i,j,0);
  double res = 0.0;

  for( int jl=-1 ; jl<=1 ; jl++ ) {
  for( int il=-1 ; il<=1 ; il++ ) {
    res += cx[il]*cy[jl] * f(i+il, j+jl, 0);
    //res += cx[il]*cy[jl] * f(ind + il + (jl*iy) );
  }}
  return res;
};


// 3D specialization
template<>
double pic::QuadraticInterpolator<3>::compute( 
        double *cx, 
        double *cy, 
        double *cz, 
        const toolbox::Mesh<float, 3>& f, 
        //const size_t ind,
        const size_t /*iy*/, 
        const size_t /*iz*/,
        int i, int j, int k)
{
  //const size_t ind = f.indx(i,j,k);
  double res = 0.0;

  for( int kl=-1 ; kl<=1 ; kl++ ){
  for( int jl=-1 ; jl<=1 ; jl++ ){
  for( int il=-1 ; il<=1 ; il++ ){
    res += cx[il]*cy[jl]*cz[kl] * f(i+il, j+jl, k+kl);
    //res += cx[il]*cy[jl]*cz[kl] * f(ind + il + (jl*iy) + (kl*iz) );
  }}}
  return res;
};



template<size_t D>
void pic::QuadraticInterpolator<D>::solve(
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
      double xpn = D >= 1 ? con.loc(0,n) - mins[0] : con.loc(0,n);
      double ypn = D >= 2 ? con.loc(1,n) - mins[1] : con.loc(1,n);
      double zpn = D >= 3 ? con.loc(2,n) - mins[2] : con.loc(2,n);

      //double u = con.vel(0,n);
      //double v = con.vel(1,n);
      //double w = con.vel(2,n);

      // particle location in the primary and dual 1/2-shifted grid

      // default 2nd order scheme
      int ip = round(xpn);
      int jp = round(ypn);
      int kp = round(zpn);
      int id = round(xpn-0.5);
      int jd = round(ypn-0.5);
      int kd = round(zpn-0.5);

      // Sokolov version; NOTE: no effect for 2nd order; maybe too narrow stencil
      //int ip = round(xpn);
      //int jp = round(ypn);
      //int kp = round(zpn);
      //int id = floor(xpn-0.5);
      //int jd = floor(ypn-0.5);
      //int kd = floor(zpn-0.5);

      //--------------------------------------------------
      // coefficients on both prime and dual (staggered +0.5) grids
      double cxd[3] = {0.0}, 
             cxp[3] = {0.0}, 
             cyd[3] = {0.0}, 
             cyp[3] = {0.0}, 
             czd[3] = {0.0}, 
             czp[3] = {0.0};
      
      // \Delta x from primary and staggered grid points
      double dxp = xpn-ip;
      double dyp = ypn-jp;
      double dzp = zpn-kp;
      double dxd = xpn-id-0.5;
      double dyd = ypn-jd-0.5;
      double dzd = zpn-kd-0.5;

      //--------------------------------------------------
      // Lorentz contract lenghts
      //double gam = sqrt(1.0 + u*u + v*v + w*w); // \gamma
      //double betax2 = pow(u/gam, 2);            // v_x^2
      //double betay2 = pow(v/gam, 2);            // v_y^2
      //double betaz2 = pow(w/gam, 2);            // v_z^2
      //double beta2  = betax2 + betay2 + betaz2; // |v|^2

      //dxp = dxp/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //dxd = dxd/(1.0 + (gam-1.0)*betax2/(beta2 + EPS));
      //dyp = dyp/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dyd = dyd/(1.0 + (gam-1.0)*betay2/(beta2 + EPS));
      //dzp = dzp/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));
      //dzd = dzd/(1.0 + (gam-1.0)*betaz2/(beta2 + EPS));

      //--------------------------------------------------
      // compute shape function weights
        
      // Default scheme
      if(D >= 1) W2nd(dxp, &cxp[0] );
      if(D >= 2) W2nd(dyp, &cyp[0] );
      if(D >= 3) W2nd(dzp, &czp[0] );
      if(D >= 1) W2nd(dxd, &cxd[0] );
      if(D >= 2) W2nd(dyd, &cyd[0] );
      if(D >= 3) W2nd(dzd, &czd[0] );

      //ver2: Eneryg conserving Sokolov alternating shape function scheme
      //W1st(dxd, &cxd[0] );
      //W1st(dyd, &cyd[0] );
      //W1st(dzd, &czd[0] );
      //W2nd(dxp, &cxp[0] );
      //W2nd(dyp, &cyp[0] );
      //W2nd(dzp, &czp[0] );

      //bool debug = false;
      //for(int iii=0; iii<3; iii++){
      //    if(cxp[iii] < 0.0)  debug = true;
      //    if(cxd[iii] < 0.0)  debug = true;

      //    if(cyp[iii] < 0.0)  debug = true;
      //    if(cyd[iii] < 0.0)  debug = true;

      //    if(czp[iii] < 0.0)  debug = true;
      //    if(czd[iii] < 0.0)  debug = true;

      //    if(cxp[iii] > 1.0)  debug = true;
      //    if(cxd[iii] > 1.0)  debug = true;
      //    if(cyp[iii] > 1.0)  debug = true;
      //    if(cyd[iii] > 1.0)  debug = true;
      //    if(czp[iii] > 1.0)  debug = true;
      //    if(czd[iii] > 1.0)  debug = true;
      //}

      //if(debug){
      //  std::cout 
      //    << "\n\ninterp xyz: " << "(" << xpn << "," << ypn << "," << zpn << ")"
      //    << " ip: "  << "(" << ip << "," << jp << "," << kp << ")"
      //    << " id: "  << "(" << id << "," << jd << "," << kd << ")"
      //    << " dxp: " << "(" << xpn - ip << "," << ypn - jp << "," << zpn - kp << ")"
      //    << " dxd: " << "(" << xpn - id + 0.5f << "," << ypn - jd + 0.5f << "," << zpn - kd + 0.5f << ")"
      //    << " W1x: " << "(" << cxp[0] << "," << cxp[1] << "," << cxp[2] << /* "," << cxp[3] << */ ")"
      //    << " W2x: " << "(" << cxd[0] << "," << cxd[1] << "," << cxd[2] << /* "," << cxd[3] << */ ")"
      //    << " W1y: " << "(" << cyp[0] << "," << cyp[1] << "," << cyp[2] << /* "," << cyp[3] << */ ")"
      //    << " W2y: " << "(" << cyd[0] << "," << cyd[1] << "," << cyd[2] << /* "," << cyd[3] << */ ")"
      //    << " W1z: " << "(" << czp[0] << "," << czp[1] << "," << czp[2] << /* "," << czp[3] << */ ")"
      //    << " W2z: " << "(" << czd[0] << "," << czd[1] << "," << czd[2] << /* "," << czd[3] << */ ")"
      //    << "\n";
      //}

      //--------------------------------------------------
      // integrate over shape function, i.e. interpolate
      // NOTE: we shift coefficient array pointers +1 to take into account -1,0,+1 access pattern
      con.ex(n) = compute( &cxd[1], &cyp[1], &czp[1], gs.ex, iy,iz,  id,jp,kp); // Ex(d,p,p)
      con.ey(n) = compute( &cxp[1], &cyd[1], &czp[1], gs.ey, iy,iz,  ip,jd,kp); // Ey(p,d,p)
      con.ez(n) = compute( &cxp[1], &cyp[1], &czd[1], gs.ez, iy,iz,  ip,jp,kd); // Ez(p,p,d)
      con.bx(n) = compute( &cxp[1], &cyd[1], &czd[1], gs.bx, iy,iz,  ip,jd,kd); // Bx(p,d,d)
      con.by(n) = compute( &cxd[1], &cyp[1], &czd[1], gs.by, iy,iz,  id,jp,kd); // By(d,p,d)
      con.bz(n) = compute( &cxd[1], &cyd[1], &czp[1], gs.bz, iy,iz,  id,jd,kp); // Bz(d,d,p)

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
template class pic::QuadraticInterpolator<1>; // 1D3V
template class pic::QuadraticInterpolator<2>; // 2D3V
template class pic::QuadraticInterpolator<3>; // 3D3V
