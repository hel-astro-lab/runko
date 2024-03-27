#include <iostream>
#include <cmath>
#include <cassert>

#include "qed/interactions/synchrotron.h"
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
  using toolbox::dot;
  using toolbox::cross;
  using toolbox::unit_cross;
  using toolbox::sum;
  using toolbox::inv;
  using toolbox::bkn_plaw;
  using toolbox::find_sorted_nearest_algo2; // binary search for sorted arrays


tuple<float, float> Synchrotron::get_minmax_ene( string t1, string t2, double ene)
{
  // only gam >1.5 emits
  return {3.0f, INF};
}


float Synchrotron::comp_chi(
    float ux1, float uy1, float uz1,
    float ex,  float ey,  float ez,
    float bx,  float by,  float bz)
{

  Vec3<float> zv; // particle 3-mom
  zv.set(ux1, uy1, uz1); 

  float gam = sqrt(1.0 + dot(zv,zv)); // prtcl gamma

  // --------------------------------------------------
  // calculate the electron quantum parameter

  Vec4<float> z(gam, ux1, uy1, uz1);  // particle four momentum

  // electromagnetic Maxwell tensor
  Vec4<float> F1( 0,  ex, ey, ez );
  Vec4<float> F2(-ex,  0,-bz, by );
  Vec4<float> F3(-ey, bz,  0, bx );
  Vec4<float> F4(-ez,-by, bx, 0  );
  Mat4<float> F(  F1, F2, F3, F4 );

  auto Fdotz = dot(F, z); // p_mu F^\mu\nu
  float chi_fpart = norm(Fdotz); //|p F| 

  //--------------------------------------------------
  // ver2
  //Vec3<float> B(bx,by,bz); // magnetic field
  //Vec3<float> E(ex,ey,ez); // electric field
  //auto beta = zv/gam; 

  //auto VxB = cross(beta, B);
  //auto E_VxB = E + VxB;
  //float chi_fpart2 = gam*sqrt( abs( pow(dot(beta, E), 2) - dot(E_VxB, E_VxB) ) );

  //std::cout << "Chi_fpart: v1:" << chi_fpart << " v2:" << chi_fpart2 << "\n";
  //std::cout << "z " << zv << "\n";
  //std::cout << "be " << beta << "\n";
  //std::cout << "B " << B << "\n";
  //std::cout << "E " << E << "\n";
  //std::cout << "vxb " << VxB << "\n";
  //std::cout << "e_vxb " << E_VxB << "\n";

  //--------------------------------------------------
  
  chi_e = chi_fpart/B_QED; // dimensionless quantum parameter \chi_pm = x
                           // NOTE we also update the internal storage for interact() function here

  return chi_e;
}


float Synchrotron::comp_optical_depth(
    string t1, 
    float ux1, float uy1, float uz1,
    float ex,  float ey,  float ez,
    float bx,  float by,  float bz)
{

  float x = comp_chi( ux1,uy1,uz1, ex,ey,ez, bx,by,bz); // dimensionless quantum parameter \chi_pm = x

  float gam = sqrt( 1.0 + ux1*ux1 + uy1*uy1 + uz1*uz1 );

  // calculate
  // K(\chi_\pm) = \int_0^\chi_\pm d^2 N_phot / dt d\chi d\chi = \int_0^\chi_pm  S(\chi_\pm, \chi_x)/\chi d\chi
  // where S is known as the Ritus synchrotron emissivity function
                                 
  // the integral K has the following asymptotic form (at x << 1 and x >> 1)
  float K_asy = 5.236f*( exp(-x) + (1.0f - exp(-x))*pow(x, -0.3333333f) ); // asymptotics

  // in addition, at x ~ 1 it has the following approximative form (found by fitting)
  // the expression is accurate to < 1%
  float K_app = 3.81108484f / (log(x/1.5746f) - ((x + x) - (-0.11058f / x)));

  //--------------------------------------------------
  // normalization of the K integral
  float prefac_int = 3.0*sqrt(3.0)*x/(4.0*PI*gam);

  // classical radiated power
  float prefac_syn = 2.0*alphaf/(3.0*lamC); 

  // total time change in the synchrotron optical depth is 
  return prefac_syn*prefac_int*(K_asy + K_app);
}


tuple<float, float> Synchrotron::accumulate(
    string t3, float e3, 
    string t4, float e4)
{

  //if( (t1 == "e-" || t1 == "e+") && e1 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl
  //if( (t2 == "e-" || t2 == "e+") && e2 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl

  float f1 = 1.0f;
  float f2 = 1.0f;

  // 1/x^a suppression below threshold
  float thr = 1.0f;
  if( e4 < thr ) f2 = pow(e4/thr, -2.0);

  return {f1,f2};
}
  


void Synchrotron::interact(
  string& t1, float& ux1, float& uy1, float& uz1,
  string& t2, float& ux2, float& uy2, float& uz2) 
{

  //--------------------------------------------------
  Vec3<float> zv, xv;
  zv.set(ux1, uy1, uz1);  // 4-velocity 3-vector
  float gam = sqrt(1.0 + dot(zv,zv)); // prtcl gamma
                                          
  //--------------------------------------------------
  // draw photon chi when electron chi is known

  // get minimum chi_x so that the photon energy is <1e-4; 
  // this sets the y-axis of XI as logspace(chi,xmin, chi_e, 33)
  float logchix_min = log10( bkn_plaw(chi_e, 3.0e-10, 1.0, 2.0, 1.0, 1.45) ); // the fitting function is accurate to <3%

  //--------------------------------------------------

  const int dim0 = 32;
  const int dim1 = 33;

  int i,j; 
  float dx, dy, logchix, chi_x;
  float XI_int[dim1]; // +1 element to ensure that the last value is always 1.0

  // closest index on the CHIE grid (x-axis of XI table)
  i = find_sorted_nearest_algo2(CHIE, chi_e, dim0);
  float rnd = rand(); // draw a random nuber

  //std::cout << "\n";
  //std::cout << "chi_e: " << chi_e << " chi_xmin: " << logchix_min << "\n";

  //--------------------------------------------------  
  if( i == 0 ) { // chi_e < chi_e,min 
    //std::cout << " XI smaller \n";

    // closest value along the y-axis (=CHIPH) 
    j = find_sorted_nearest_algo2(XI[0], rnd, dim1);
    dy = log10( XI[0][j] ) - log10( rnd );
    dy /= abs( log10(XI[0][j]) - log10(XI[0][j-1]) ); // normalize

  //--------------------------------------------------  
  } else if (i == dim0) { // chi_e > chi_max
    //std::cout << " XI larger \n";

    // closest value along the y-axis (=CHIPH) 
    j = find_sorted_nearest_algo2(XI[dim0-1], rnd, dim1);
    dy = log10(XI[dim0-1][j]) - log10(rnd);
    dy /= abs( log10(XI[dim0-1][j]) - log10(XI[dim0-1][j-1]) ); // normalize

  //--------------------------------------------------  
  } else { // chi_e,min < chi_e < chi_e,max; interpolate
    //std::cout << " interpolating XI \n";

    dx = log10(CHIE[i]) - log10(chi_e);  // log difference 
    dx /= abs( log10(CHIE[i]) - log10(CHIE[i-1]) ); // normalize

    //std::cout << " chi_e:" << chi_e << "\n";
    //std::cout << " i:" << i << " CHIE[i]" << CHIE[i] << " dx:" << dx << "\n";

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
  float dy_logchix = ( log10(chi_e) - logchix_min )/(dim1 - 1);

  if( j == 0 ){
    dy = 0.0;
    logchix = logchix_min; 

  //} else if (j==dim1) {
    // NOTE this branch is never reached because last point is always 1.0

  } else { // interpolate in y-dim

    // logarithmic value from chi_x,min to chi_e
    logchix = logchix_min + dy_logchix*(j-1 + dy); 
  }

  //std::cout << " rnd:" << rnd << "\n";
  //std::cout << " j:" << j << " XI[j]" << XI_int[j] << " CHIE[j]:" << CHIE[j] << " dy:" << dy << "\n";

  chi_x = pow(10.0f, logchix);

  //std::cout << " CHIX:" << chi_x << " CHIE: " << chi_e << "\n";

  //--------------------------------------------------
  // construct photon properties 
    
  // photon energy
  float x = (gam-1.0f)*chi_x/chi_e; // kinetic energy  
  //float x = norm(zv)*chi_x/chi_e; // momentum \beta\gamma 

  //std::cout << " xsyn:" << x << " game" << gam << " xsyn/(gam-1):" << x/(gam-1.0f) << "\n";

  // assert that photon energy is reasonable
  //assert( x < gam - 1.0f );

  float z0 = norm(zv); // length of the electron velocity vector

  // photon to the direction of the electron 
  t2 = "ph";
  ux2 = x*zv(0)/z0;
  uy2 = x*zv(1)/z0;
  uz2 = x*zv(2)/z0;

  //--------------------------------------------------
  // udpate electron properties

  // remove photon momentum from electron
  ux1 -= ux2;
  uy1 -= uy2;
  uz1 -= uz2;

  return;
}


} // end of namespace qed

#pragma GCC diagnostic pop
