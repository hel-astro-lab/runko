#include "quadratic_2nd.h"

#include <cmath> 
#include <cassert>

#include "../../tools/iter/iter.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif


// Cloud-in-cell shape
// 1st order weight
//inline void compute_coeffs(double d, double* coeff){
//  coeff[0] = 0.0;
//  coeff[1] = 1.0-d;
//  coeff[2] = d;
//}


// Triangular shaped cloud;
// 2nd order quadratic spline
/*      -
*       |  3/4 - x^2                  if |x|<1/2
* W(x)=<   1/2*(3/2 - |x|)^2          if 1/2<=|x|<3/2
*       |  0                          otherwise
*       -
*/
//9/8 (1- 2/3 x)^2
//(3/4- 1/2 x)^2
//1/4*(3/2-x)**2
//
  // TODO: c0 and c2 appear wrong; should be 
  //   1/8 (3-2x)^2 
  // = 0.5*( 1.5-x )^2  
  // = 0.5*( x^2 - 3x + 2.25)
inline void compute_coeffs(double d, double* coeff){

  //NOTE: d at wings includes +-1 so W2_mp1 -> 0.5*(3/2 +- d)^2
  coeff[0] = 0.50*(0.5 - d)*(0.5 - d); //W2_im1 
  coeff[1] = 0.75 - d*d;                 //W2_i   
  coeff[2] = 0.50*(0.5 + d)*(0.5 + d); //W2_ip1 

}


// 1st order weight for Sokolov's alternating scheme
inline void compute_lower_coeffs(double d, double* coeff){

  coeff[0] = 0.0;
  coeff[1] = 1.0-d;
  coeff[2] = d;

  // Sokolov higher
  // 0.5*(1-d)^2 ok
  // 3/4 - (0.5 - d)^2
  // d^2
}


// Cubic piewce spline  
// 3rd order  
/*
inline void compute_coeffs(float_m dx, float_m* coeff){
  float_m dx2 = dx*dx;
  float_m dx3 = dx*dx*dx;

  coeff[0] = ( 1.0f - dx3 )*1.0f/6.0f - 0.5f*( dx-dx2 );
  coeff[1] = 2.0f/3.0f - xi2 + 0.5f*dx3;
  coeff[2] = 1.0f/6.0f + 0.5f*( dx+dx2-dx3 );
  coeff[3] = dx3*1.0f/6.0f;
}
*/

// Cubic piewce spline  
// 4th order  
/*
inline void compute_coeffs(float_m dx, float_m* coeff){
  float_m dx2 = dx*dx;
  float_m dx3 = dx*dx*dx;
  float_m dx4 = dx*dx*dx*dx;

  coeff[0] = 1.0f  /384.0f - 1.0f /48.0f*dx  + 1.0f/16.0f*dx2 - 1.0f/12.0f*dx3 + 1.0f/24.0f*dx4;
  coeff[1] = 19.0f /96.0f  - 11.0f/24.0f*dx  + 1.0f/4.0f* dx2 + 1.0f/6.0f* dx3 - 1.0f/6.0f* dx4;
  coeff[2] = 115.0f/192.0f - 5.0f /8.0f *dx2 + 1.0f/4.0f* dx4;
  coeff[3] = 19.0f /96.0f  + 11.0f/24.0f*dx  + 1.0f/4.0f* dx2 - 1.0f/6.0f* dx3 - 1.0f/6.0f* dx4;
  coeff[4] = 1.0f  /384.0f + 1.0f /48.0f*dx  + 1.0f/16.0f*dx2 + 1.0f/12.0f*dx3 + 1.0f/24.0f*dx4;
}
*/

// Dummy template default 
template<std::size_t D>
inline double compute( double* /*cx*/, double* /*cy*/, double* /*cz*/, 
        const toolbox::Mesh<float_m, 3>& /*f*/, 
        const size_t /*iy*/, const size_t /*iz*/,
        int /*i*/, int /*j*/, int /*k*/,
        std::array<int, D>& /*dummy template specialization variable*/
        )
{ return 0.0; }


// 2D specialization
template<>
inline double compute( 
        double* cx, 
        double* cy, 
        double* /*cz*/, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t iy, 
        const size_t /*iz*/,
        int i, int j, int /*k*/,
        std::array<int, 2>& /*dummy template specialization variable*/
        )
{
  const size_t ind = f.indx(i,j,0);
  float_p res = 0.0f;

  for( int jl=-1 ; jl<=1 ; jl++ ) {
  for( int il=-1 ; il<=1 ; il++ ) {
    res += cx[il]*cy[jl] * f(ind + il + (jl*iy) );
  }}
  return res;
};


// 3D specialization
template<>
inline double compute( 
        double *cx, 
        double *cy, 
        double *cz, 
        const toolbox::Mesh<float_m, 3>& f, 
        //const size_t ind,
        const size_t iy, 
        const size_t iz,
        int i, int j, int k,
        std::array<int, 3>& /*dummy template specialization variable*/
        )
{
  const size_t ind = f.indx(i,j,k);
  float_p res = 0.0f;

  for( int kl=-1 ; kl<=1 ; kl++ ) {
  for( int jl=-1 ; jl<=1 ; jl++ ) {
  for( int il=-1 ; il<=1 ; il++ ) {
    //res += cx[il]*cy[jl]*cz[kl] * f(idx + il, idy + jl, idz + kl);
    res += cx[il]*cy[jl]*cz[kl] * f(ind + il + (jl*iy) + (kl*iz) );
  }}}
  return res;
};



template<size_t D, size_t V>
void pic::QuadraticInterpolator<D,V>::solve(
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

      //compute_coeffs( xpn - ip,        &cxp[0] );
      //compute_coeffs( ypn - jp,        &cyp[0] );
      //compute_coeffs( zpn - kp,        &czp[0] );
      //compute_coeffs( xpn - id + 0.5f, &cxd[0] );
      //compute_coeffs( ypn - jd + 0.5f, &cyd[0] );
      //compute_coeffs( zpn - kd + 0.5f, &czd[0] );
        
      
      // \Delta x from primary and staggered grid points
      double dxp = xpn-ip;
      double dyp = ypn-jp;
      double dzp = zpn-kp;

      double dxd = xpn-id-0.5;
      double dyd = ypn-jd-0.5;
      double dzd = zpn-kd-0.5;
        
      //double gam = sqrt(1.0 + u*u + v*v + w*w);

      //DONE ver2: DONE: sokolov
      if(D >= 1) compute_lower_coeffs(dxp, &cxp[0] );
      if(D >= 2) compute_lower_coeffs(dyp, &cyp[0] );
      if(D >= 3) compute_lower_coeffs(dzp, &czp[0] );
      if(D >= 1) compute_coeffs(      dxd, &cxd[0] );
      if(D >= 2) compute_coeffs(      dyd, &cyd[0] );
      if(D >= 3) compute_coeffs(      dzd, &czd[0] );

      // default
      //if(D >= 1) compute_coeffs(dxp, &cxp[0] );
      //if(D >= 2) compute_coeffs(dyp, &cyp[0] );
      //if(D >= 3) compute_coeffs(dzp, &czp[0] );
      //if(D >= 1) compute_coeffs(dxd, &cxd[0] );
      //if(D >= 2) compute_coeffs(dyd, &cyd[0] );
      //if(D >= 3) compute_coeffs(dzd, &czd[0] );


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
      con.ex(n) = compute<D>( &cxd[1], &cyp[1], &czp[1], yee.ex, iy,iz,  id,jp,kp, dummy); // Ex(d,p,p)
      con.ey(n) = compute<D>( &cxp[1], &cyd[1], &czp[1], yee.ey, iy,iz,  ip,jd,kp, dummy); // Ey(p,d,p)
      con.ez(n) = compute<D>( &cxp[1], &cyp[1], &czd[1], yee.ez, iy,iz,  ip,jp,kd, dummy); // Ez(p,p,d)
      con.bx(n) = compute<D>( &cxp[1], &cyd[1], &czd[1], yee.bx, iy,iz,  ip,jd,kd, dummy); // Bx(p,d,d)
      con.by(n) = compute<D>( &cxd[1], &cyp[1], &czd[1], yee.by, iy,iz,  id,jp,kd, dummy); // By(d,p,d)
      con.bz(n) = compute<D>( &cxd[1], &cyd[1], &czp[1], yee.bz, iy,iz,  id,jd,kp, dummy); // Bz(d,d,p)

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

template class pic::QuadraticInterpolator<2,3>; // 2D3V
template class pic::QuadraticInterpolator<3,3>; // 3D3V
