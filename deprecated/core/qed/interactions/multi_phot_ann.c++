#include <iostream>
#include <cmath>
#include <cassert>

#include "core/qed/interactions/multi_phot_ann.h"
#include "tools/vector.h"
#include "tools/sample_arrays.h"
#include "tools/bkn_plaw.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace qed {

  using std::string;
  using std::tuple;
  using std::sqrt;

  using toolbox::Vec3;
  using toolbox::Mat3;
  using toolbox::Vec4;
  using toolbox::Mat4;
  using toolbox::norm;
  using toolbox::norm_minkowski;
  using toolbox::dot;
  using toolbox::cross;
  using toolbox::unit_cross;
  using toolbox::sum;
  using toolbox::inv;
  using toolbox::bkn_plaw;
  using toolbox::find_sorted_nearest_algo2; // binary search for sorted arrays



tuple<float, float> MultiPhotAnn::get_minmax_ene( string /*t1*/, string /*t2*/, double /*ene*/)
{
  // only x>2 can participate
  return {2.0f, INF};
}


float MultiPhotAnn::comp_chi(
    float ux1, float uy1, float uz1,
    float ex,  float ey,  float ez,
    float bx,  float by,  float bz)
{

  Vec3<float> xv; // particle 3-mom
  xv.set(ux1, uy1, uz1); 
  float x0 = norm(xv); // photon energy

  // --------------------------------------------------
  // calculate the photon quantum parameter

  Vec4<float> z(x0, ux1, uy1, uz1);  // particle four momentum

  // contravariant electromagnetic Maxwell tensor F^\mu\nu
  Vec4<float> F1( 0, -ex,-ey,-ez );
  Vec4<float> F2( ex,  0,-bz, by );
  Vec4<float> F3( ey, bz,  0,-bx );
  Vec4<float> F4( ez,-by, bx, 0  );
  Mat4<float> F(  F1, F2, F3, F4 );

  // covariant electromagnetic Maxwell tensor F_\mu\nu
  // metric signature of +---
  //Vec4<float> F1( 0,  ex, ey, ez );
  //Vec4<float> F2(-ex,  0,-bz, by );
  //Vec4<float> F3(-ey, bz,  0,-bx );
  //Vec4<float> F4(-ez,-by, bx, 0  );
  //Mat4<float> F(  F1, F2, F3, F4 );

  auto Fdotz = dot(F, z); // p_mu F^\mu\nu
  float chi_fpart = norm_minkowski(Fdotz); //|p F| 


  ////--------------------------------------------------
  /// NOTE for implementation comparisons, see synchrotorn.c++

  //--------------------------------------------------
  chi_x = chi_fpart/B_QED; // dimensionless quantum parameter \chi_x 
                           // NOTE we also update the internal storage for interact() function here

  return chi_x;
}


float MultiPhotAnn::comp_optical_depth(
    string /*t1*/, 
    float ux1, float uy1, float uz1,
    float ex,  float ey,  float ez,
    float bx,  float by,  float bz)
{
    
  float xene = sqrt( ux1*ux1 + uy1*uy1 + uz1*uz1 ); // photon energy

  float x = comp_chi( ux1,uy1,uz1, ex,ey,ez, bx,by,bz); // dimensionless quantum parameter \chi_x

  // asymptotic behavior of the T function
  float T_asy = 1.68f*exp(-2.8f/x)*pow(x,1.7f);

  // TODO real analytical asymptotic behavior is
  // 1.249 exp(-8/3/x)*x^2 and 2.0678*x^(5/3)
  // need to fix this and re-fit

  //--------------------------------------------------
  // empirical fitting function to correct the behavior at x ~ 1
  // fit is accurate to <0.2% for x > 1

  // log-gaussian
  const float a1 = 0.275;
  const float mu1 = -0.89846;
  const float sig1 = 7.6766;
  float T_corr = 1.0f - a1*exp(-pow(log(x)-mu1, 2.0f)/sig1);

  // normalization of the T integral
  float prefac_int = 1.0/(PI*sqrt(3.0)*x*xene);
    
  // classical BBW power
  float prefac_bbw = alphaf/lamC; 
    
  return prefac_bbw*prefac_int*T_asy*T_corr;
}



void MultiPhotAnn::interact(
  string& t1, float& ux1, float& uy1, float& uz1,
  string& t2, float& ux2, float& uy2, float& uz2) 
{

  //--------------------------------------------------
  Vec3<float> xv, zm, zp;
  xv.set(ux1, uy1, uz1);  // 4-velocity 3-vector
  float x0 = norm(xv); // photon gamma
                                          
  //--------------------------------------------------
  // draw electron chi when photon chi is known

  //std::cout << " \n";

  // get minimum chi_e so that the photon energy is <1e-4; 
  // this sets the y-axis of XI as logspace(chi,emin, chi_x, 33)
  // the fitting function is accurate to <0.2%
  float logchie_min = log10( bkn_plaw(chi_x, 0.01876, 0.103, 0.931, 0.0475, 1.08) );
                                                                                
  //--------------------------------------------------
  const int dim0 = 32;
  const int dim1 = 33;

  int i,j; 
  float dx=0, dy, logchie, chi_ep;
  float XI_int[dim1]; // +1 element to ensure that the last value is always 1.0

  // closest index on the CHIE grid (x-axis of XI table)
  i = find_sorted_nearest_algo2(CHIX, chi_x, dim0);
  float rnd = rand(); // draw a random number

  //--------------------------------------------------  
  if( i == 0 ) { // chi_x < chi_x,min 
    //std::cout << " XI smaller \n";

    // closest value along the y-axis (\propto CHIE) 
    j = find_sorted_nearest_algo2(XI[0], rnd, dim1);
    dy = log10( XI[0][j] ) - log10( rnd );
    dy /= abs( log10(XI[0][j]) - log10(XI[0][j-1]) ); // normalize
                                                      //
  //--------------------------------------------------  
  } else if (i == dim0) { // chi_x > chi_max
    //std::cout << " XI larger \n";

    // closest value along the y-axis (\propto CHIX) 
    j = find_sorted_nearest_algo2(XI[dim0-1], rnd, dim1);
    dy = log10(XI[dim0-1][j]) - log10(rnd);
    dy /= abs( log10(XI[dim0-1][j]) - log10(XI[dim0-1][j-1]) ); // normalize
                                                                  
  //--------------------------------------------------  
  } else { // chi_x,min < chi_x < chi_x,max; interpolate
    //std::cout << " interpolating XI \n";

    dx = log10(CHIX[i]) - log10(chi_x);  // log difference 
    dx /= abs( log10(CHIX[i]) - log10(CHIX[i-1]) ); // normalize

    //std::cout << " chi_x:" << chi_x << "\n";
    //std::cout << " i:" << i << " CHIX[i]" << CHIX[i] << " dx:" << dx << "\n";

    assert(dx >= 0.0f);
    assert(dx <= 1.0f);

    // interpolate in x-dim; done in log10-log10 space; 
    // TODO it would a tiny bit be faster to do this in e-base
    for(size_t n=0; n<dim1; n++) XI_int[n] = pow(10.0f, (1.0f-dx)*log10(XI[i-1][n]) + dx*log10(XI[i][n]));

    // closest value along the y-axis (=CHIPH) 
    j = find_sorted_nearest_algo2(XI_int, rnd, dim1);

    dy = log10(XI_int[j]) - log10(rnd); // this can overflow; however, in that case it is not used 
    dy /= abs( log10(XI_int[j]) - log10(XI_int[j-1]) ); // normalize
  }

  //--------------------------------------------------
  // we have found the correct x-slice; now find y-value
  float dy_logchie = ( log10(chi_x) - logchie_min )/(dim1 - 1);


  if( j == 0 ){
    dy = 0.0;
    logchie = logchie_min; 

  } else { // interpolate in y-dim

    // logarithmic value from chi_e,min to chi_x
    logchie = logchie_min + dy_logchie*(j-1 + dy); 
  }

  //std::cout << " rnd:" << rnd << "\n";
  //std::cout << " j:" << j << " XI[j]" << XI_int[j] << " CHIX[j]:" << CHIX[j] << " dy:" << dy << "\n";
    
  chi_ep = pow(10.0f, logchie);

  //--------------------------------------------------
  // xi table is symmetric around maximum value; need to randomize which particle gets more energy

  float chi_e, chi_p; // electron and positron quantum numbers
  if( rand() < 0.5 ) {
    chi_e = chi_ep;
    chi_p = chi_x - chi_e;
  } else {
    chi_p = chi_ep;
    chi_e = chi_x - chi_p;
  }

  // cap small values; the substraction can sometimes result in negative values when chi_e ~ chi_x
  chi_e = chi_e < 0.0 ? 1.0e-4 : chi_e;
  chi_p = chi_p < 0.0 ? 1.0e-4 : chi_p;

  //--------------------------------------------------
  // particles are emtited to the direction of the photon
  // TODO: better way would be to go into the photons frame, draw the angles, create leptons, and de-boost back

  float inv_cx = ( x0 - 2.0 )/chi_x; // available energy / chi_x
    
  t1 = "e-";
  float pe = sqrt( pow( 1.0 + chi_e*inv_cx, 2 ) - 1.0 );
  ux1 = pe*xv(0)/x0; 
  uy1 = pe*xv(1)/x0; 
  uz1 = pe*xv(2)/x0; 

  t2 = "e+";
  float pp = sqrt( pow( 1.0 + chi_p*inv_cx, 2 ) - 1.0 );
  ux2 = pp*xv(0)/x0; 
  uy2 = pp*xv(1)/x0; 
  uz2 = pp*xv(2)/x0; 
                                                                                

  //--------------------------------------------------
#ifdef DEBUG 
  
  // TODO check momentum conservation

  // check for NaNs
  bool ts[6];
  ts[0] = std::isnan(ux1);
  ts[1] = std::isnan(uy1);
  ts[2] = std::isnan(uz1);
  ts[3] = std::isnan(ux2);
  ts[4] = std::isnan(uy2);
  ts[5] = std::isnan(uz2);

  if(ts[0] ||
     ts[1] ||
     ts[2] ||
     ts[3] ||
     ts[4] ||
     ts[5] ) { 

    std::cerr << "ERROR MULTI-PHOT-ANN:" << std::endl;
    for(size_t i = 0; i < 5; i++) { std::cerr << i << " " << ts[i] << std::endl; }

    std::cerr << "logchie:" << logchie << " chi_ep:" << chi_ep <<  " chi_e:" <<  chi_e << " chi_p:" << chi_p << std::endl;
    std::cerr << "xv: " <<  xv << std::endl;
    std::cerr << "x0: " <<  x0 << " inv_cx:" << inv_cx << std::endl;
    std::cerr << "pp: " <<  pp << " pe:" << pe << std::endl;

    std::cerr << "ux1: " <<  ux1 << " uy1:" << uy1 << " uz1:" << uz1 << std::endl;
    std::cerr << "ux2: " <<  ux2 << " uy2:" << uy2 << " uz2:" << uz2 << std::endl;

    std::cerr << "i: " <<  i << " dx:" << dx << std::endl;
    std::cerr << "j: " <<  j << " dy:" << dy << std::endl;
  }

#endif

  return;
}


} // end of namespace qed
    
#pragma GCC diagnostic pop
