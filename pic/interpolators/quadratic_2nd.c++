#include "quadratic_2nd.h"

#include <cmath> 
#include <cassert>

#include "../shapes.h"
#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif



// 2D specialization
template<>
double pic::QuadraticInterpolator<2>::compute( 
        double* cx, 
        double* cy, 
        double* /*cz*/, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t iy, 
        const size_t /*iz*/,
        int i, int j, int /*k*/)
{
  const size_t ind = f.indx(i,j,0);
  double res = 0.0;

  for( int jl=-1 ; jl<=1 ; jl++ ) {
  for( int il=-1 ; il<=1 ; il++ ) {
    res += cx[il]*cy[jl] * f(ind + il + (jl*iy) );
  }}
  return res;
};


// 3D specialization
template<>
double pic::QuadraticInterpolator<3>::compute( 
        double *cx, 
        double *cy, 
        double *cz, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t iy, 
        const size_t iz,
        int i, int j, int k)
{
  const size_t ind = f.indx(i,j,k);
  double res = 0.0;

  for( int kl=-1 ; kl<=1 ; kl++ ) {
  for( int jl=-1 ; jl<=1 ; jl++ ) {
  for( int il=-1 ; il<=1 ; il++ ) {
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

  //NOTE: we need this dummy to enforce template D specialization later on
  std::array<int, D> dummy; 

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
    //            fields::YeeLattice& yee,
    //            pic::ParticleContainer<D>& con){
    for(size_t n=0; n<con.size(); n++) {
      //--------------------------------------------------
      // indices & locs to grid
      double xpn = con.loc(0, n) - mins[0];
      double ypn = con.loc(1, n) - mins[1];
      double zpn = con.loc(2, n) - mins[2];

      double u = con.vel(0,n);
      double v = con.vel(1,n);
      double w = con.vel(2,n);


      // particle location in the primary and dual 1/2-shifted grid
      // TODO: round or floor (=corrected int() for neg numbers)
      // TODO: if(D >= 1) switches
      // TODO: stagger up or down?
      
      // TODO ver
      //int ip = round(xpn);
      //int jp = round(ypn);
      //int kp = round(zpn);
      //int id = round(xpn + 0.5f);
      //int jd = round(ypn + 0.5f);
      //int kd = round(zpn + 0.5f);

      // TODO ver2
      int ip = floor(xpn)-1.0f;
      int jp = floor(ypn)-1.0f;
      int kp = floor(zpn)-1.0f;
      int id = floor(xpn);
      int jd = floor(ypn);
      int kd = floor(zpn);

      //TODO ver3
      //int ip = floor(xpn);
      //int jp = floor(ypn);
      //int kp = floor(zpn);
      //int id = floor(xpn)+1.0f;
      //int jd = floor(ypn)+1.0f;
      //int kd = floor(zpn)+1.0f;


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
        
      //ver2: Eneryg conserving Sokolov alternating shape function scheme
      if(D >= 1) W1st(dxp, &cxp[0] );
      if(D >= 2) W1st(dyp, &cyp[0] );
      if(D >= 3) W1st(dzp, &czp[0] );

      if(D >= 1) W2nd(dxd, &cxd[0] );
      if(D >= 2) W2nd(dyd, &cyd[0] );
      if(D >= 3) W2nd(dzd, &czd[0] );

      // default scheme
      //if(D >= 1) W2nd(dxp, &cxp[0] );
      //if(D >= 2) W2nd(dyp, &cyp[0] );
      //if(D >= 3) W2nd(dzp, &czp[0] );
      //if(D >= 1) W2nd(dxd, &cxd[0] );
      //if(D >= 2) W2nd(dyd, &cyd[0] );
      //if(D >= 3) W2nd(dzd, &czd[0] );


      //std::cout 
      //    << " xyz: " << "(" << xpn << "," << ypn << "," << zpn << ")"
      //    << " ip: "  << "(" << ip << "," << jp << "," << kp << ")"
      //    << " id: "  << "(" << id << "," << jd << "," << kd << ")"
      //    << " dxp: " << "(" << xpn - ip << "," << ypn - jp << "," << zpn - kp << ")"
      //    << " dxd: " << "(" << xpn - id + 0.5f << "," << ypn - jd + 0.5f << "," << zpn - kd + 0.5f << ")"
      //    << "\n";

      //--------------------------------------------------
      // integrate over shape function, i.e. interpolate
      // NOTE: we shift coefficient array pointers +1 to take into account -1,0,+1 access pattern
      con.ex(n) = compute( &cxd[1], &cyp[1], &czp[1], yee.ex, iy,iz,  id,jp,kp); // Ex(d,p,p)
      con.ey(n) = compute( &cxp[1], &cyd[1], &czp[1], yee.ey, iy,iz,  ip,jd,kp); // Ey(p,d,p)
      con.ez(n) = compute( &cxp[1], &cyp[1], &czd[1], yee.ez, iy,iz,  ip,jp,kd); // Ez(p,p,d)
      con.bx(n) = compute( &cxp[1], &cyd[1], &czd[1], yee.bx, iy,iz,  ip,jd,kd); // Bx(p,d,d)
      con.by(n) = compute( &cxd[1], &cyp[1], &czd[1], yee.by, iy,iz,  id,jp,kd); // By(d,p,d)
      con.bz(n) = compute( &cxd[1], &cyd[1], &czp[1], yee.bz, iy,iz,  id,jd,kp); // Bz(d,d,p)

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

template class pic::QuadraticInterpolator<2>; // 2D3V
template class pic::QuadraticInterpolator<3>; // 3D3V
