#include <iostream>
#include <cmath>
#include <cassert>
#include "multi_phot_ann.h"
#include "../../tools/vector.h"
#include "../../tools/sample_arrays.h"
#include "../../tools/bkn_plaw.h"


namespace qed {

  using std::string;
  using std::tuple;
  using std::sqrt;

  using toolbox::Vec3;
  using toolbox::Mat3;
  using toolbox::Vec4;
  using toolbox::Mat4;
  using toolbox::norm;
  using toolbox::dot;
  using toolbox::cross;
  using toolbox::unit_cross;
  using toolbox::sum;
  using toolbox::inv;
  using toolbox::bkn_plaw;
  using toolbox::find_sorted_nearest_algo2; // binary search for sorted arrays



tuple<float_p, float_p> MultiPhotAnn::get_minmax_ene( string t1, string t2, double ene)
{
  // only x>2 can participate
  return {2.0f, INF};
}


float_p MultiPhotAnn::comp_chi(
    float_p ux1, float_p uy1, float_p uz1,
    float_p ex,  float_p ey,  float_p ez,
    float_p bx,  float_p by,  float_p bz)
{

  Vec3<float_p> xv; // particle 3-mom
  xv.set(ux1, uy1, uz1); 
  float_p x0 = norm(xv); // photon energy

  // --------------------------------------------------
  // calculate the photon quantum parameter

  Vec4<float_p> z(x0, ux1, uy1, uz1);  // particle four momentum

  // electromagnetic Maxwell tensor
  Vec4<float_p> F1( 0,  ex, ey, ez );
  Vec4<float_p> F2(-ex,  0,-bz, by );
  Vec4<float_p> F3(-ey, bz,  0, bx );
  Vec4<float_p> F4(-ez,-by, bx, 0  );
  Mat4<float_p> F(  F1, F2, F3, F4 );

  auto Fdotz = dot(F, z); // p_mu F^\mu\nu
  float_p chi_fpart = norm(Fdotz); //|p F| 


  //--------------------------------------------------
  chi_x = chi_fpart/B_QED; // dimensionless quantum parameter \chi_x 
                           // NOTE we also update the internal storage for interact() function here

  return chi_x;
}


float_p MultiPhotAnn::comp_optical_depth(
    string t1, 
    float_p ux1, float_p uy1, float_p uz1,
    float_p ex,  float_p ey,  float_p ez,
    float_p bx,  float_p by,  float_p bz)
{
    

  float_p x = comp_chi( ux1,uy1,uz1, ex,ey,ez, bx,by,bz); // dimensionless quantum parameter \chi_x

  // asymptotic behavior of the T function
  float_p T_asy = 1.68f*exp(-2.8f/x)*pow(x,1.7f);

  // TODO real analytical asymptotic behavior is
  // 1.249 exp(-8/3/x) and 2.0678*x^(5/3)
  // need to fix this and re-fit

  //--------------------------------------------------
  // empirical fitting function to correct the behavior at x ~ 1
  // fit is accurate to <0.2% for x > 1

  // log-gaussian
  const float_p a1 = 0.275;
  const float_p mu1 = -0.89846;
  const float_p sig1 = 7.6766;
  float_p T_corr = 1.0f - a1*exp(-pow(log(x)-mu1, 2.0f)/sig1);

  // TODO add numerical prefactor

  //float_p prefac = alpha/(lamC * PI*sqrt(3.0) );
    
  return T_asy*T_corr;
}



void MultiPhotAnn::interact(
  string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
  string& t2, float_p& ux2, float_p& uy2, float_p& uz2) 
{

  //--------------------------------------------------
  Vec3<float_p> xv, zm, zp;
  xv.set(ux1, uy1, uz1);  // 4-velocity 3-vector
  float_p x0 = norm(xv); // photon gamma
                                          
  //--------------------------------------------------
  // draw electron chi when photon chi is known

  //std::cout << " \n";

  // get minimum chi_e so that the photon energy is <1e-4; 
  // this sets the y-axis of XI as logspace(chi,emin, chi_x, 33)
  // the fitting function is accurate to <0.2%
  float_p logchie_min = log10( bkn_plaw(chi_x, 0.01876, 0.103, 0.931, 0.0475, 1.08) );
                                                                                
  //--------------------------------------------------
  const int dim0 = 32;
  const int dim1 = 33;

  int i,j; 
  float dx, dy, logchie, chi_ep;
  float XI_int[dim1]; // +1 element to ensure that the last value is always 1.0

  // closest index on the CHIE grid (x-axis of XI table)
  i = find_sorted_nearest_algo2(CHIX, chi_x, dim0);
  float rnd = rand(); // draw a random nuber

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
  float_p dy_logchie = ( log10(chi_x) - logchie_min )/(dim1 - 1);


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

  float_p chi_e, chi_p; // electron and positron quantum numbers
  if( rand() < 0.5 ) {
    chi_e = chi_ep;
    chi_p = chi_x - chi_e;
  } else {
    chi_p = chi_ep;
    chi_e = chi_x - chi_p;
  }

  //--------------------------------------------------

  //std::cout << " chi x e p" << chi_x << " " << chi_e << " " << chi_p << "\n";


  // particles are emtited to the direction of the photon
    
  // TODO: better way would be to go into the photons frame, draw the angles, create leptons, and de-boost back

  float_p inv_cx = ( x0 - 2.0 )/chi_x; // available energy / chi_x
    
  t1 = "e-";
  float_p pe = sqrt( pow( 1.0 + chi_e*inv_cx, 2 ) - 1.0 );
  ux1 = pe*xv(0)/x0; 
  uy1 = pe*xv(1)/x0; 
  uz1 = pe*xv(2)/x0; 

  t2 = "e+";
  float_p pp = sqrt( pow( 1.0 + chi_p*inv_cx, 2 ) - 1.0 );
  ux2 = pp*xv(0)/x0; 
  uy2 = pp*xv(1)/x0; 
  uz2 = pp*xv(2)/x0; 
                                                                                
  
  float_p err_ene = x0 - (sqrt(1.0 + pe*pe) + sqrt(1.0 + pp*pp));

  //std::cout << " enes x/e/p" << x0 << " " << pe << " " << pp << " err:" << err_ene << "\n";
  
  // TODO check momentum conservation


  return;
}


} // end of namespace qed
