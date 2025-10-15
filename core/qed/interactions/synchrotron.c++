#include <iostream>
#include <cmath>
#include <cassert>

#include "core/qed/interactions/synchrotron.h"
#include "tools/vector.h"
#include "tools/sample_arrays.h"
#include "tools/bkn_plaw.h"
#include "tools/signum.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace qed {

  using std::string;
  using std::tuple;
  using std::sqrt;
  using std::cbrt;

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
  using toolbox::sign;


tuple<float, float> Synchrotron::get_minmax_ene( string /*t1*/, string /*t2*/, double /*ene*/)
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

  Vec4<float> z(gam, ux1, uy1, uz1);  // particle four momentum z^mu

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

  // TODO why does this not match 

  auto Fdotz = dot(F, z); // p_mu F^\mu\nu
  float chi_fpart = norm_minkowski(Fdotz); //|p F| 

  //--------------------------------------------------
  //if(chi_fpart/B_QED > 0.04) {
  //  std::cout << " chi sy v1: " << chi_fpart 
  //            << " chie:" << chi_fpart/B_QED 
  //            << " Bq:" << B_QED 
  //            << " u:" << z 
  //            << " F.u:" << Fdotz 
  //            << " Ex/Bq" << ex/B_QED 
  //            << " By/Bq" << gam*by/B_QED 
  //            << "\n";
  //}


  // DONE manual calc (v4 with covariant F) agrees with manual version (v2)
  //      note the transpose when compiling F, hence the co <-> contra flip

  ////--------------------------------------------------
  //// ver2
  //Vec3<float> B(bx,by,bz); // magnetic field
  //Vec3<float> E(ex,ey,ez); // electric field
  //auto beta = zv/gam; 

  //auto ExB = cross(E, B);
  //auto VxB = cross(beta, B);
  //auto E_VxB = E + VxB;
  //auto v_E = dot(beta, E);
  //float chi_fpart4 = gam*sqrt( abs( v_E*v_E - dot(E_VxB, E_VxB) ) );

  //std::cout << " chi sy v4: " << chi_fpart4/B_QED << "\n";

  ////std::cout << "Chi_fpart: v1:" << chi_fpart << " v2:" << chi_fpart4 << "\n";
  //std::cout << "gam " << gam << "\n";
  //std::cout << "p_mu F^mu nu " << Fdotz << "\n";
  //std::cout << "z " << zv << "\n";
  //std::cout << "beta " << beta << "\n";
  //std::cout << "B " << B << "\n";
  //std::cout << "E " << E << "\n";
  //std::cout << "VxB " << VxB << "\n";
  //std::cout << "e_vxb " << E_VxB << "\n";
  //std::cout << "exb " << ExB << "\n";

  // TODO calculate pitch angle next

  //--------------------------------------------------
  
  chi_e = chi_fpart/B_QED; // dimensionless quantum parameter \chi_pm = x
                           // NOTE we also update the internal storage for interact() function here

  return chi_e;
}


float Synchrotron::comp_optical_depth(
    string /*t1*/, 
    float ux1, float uy1, float uz1,
    float ex,  float ey,  float ez,
    float bx,  float by,  float bz)
{

  float x = comp_chi( ux1,uy1,uz1, ex,ey,ez, bx,by,bz); // dimensionless quantum parameter \chi_pm = x

  float gam = sqrt( 1.0 + ux1*ux1 + uy1*uy1 + uz1*uz1 );

  //std::cout << "  gam:" << gam << " chi:" << x << " chi/gam:" << x/gam << std::endl;

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
  float prefac_syn = 2.0*alphaf*cvel/(3.0*lamC); 

  // total time change in the synchrotron optical depth is 
  //return prefac_syn*prefac_int*(K_asy + K_app);
  float dtau_dt = prefac_syn*prefac_int*(K_asy + K_app);


  //--------------------------------------------------
  // v2 with cleaner syntax
  //float prefac = alphaf/lamC;
  //float T = 1.442;
  //float dtau_dt2 = prefac*x*T/gam;

  float dtau_dt2 = (alphaf*cvel/lamC)*1.442*x/gam;

  //std::cout << "dtau_dt: " << dtau_dt << " " << dtau_dt2 << std::endl;

  return dtau_dt2;
}


tuple<float, float> Synchrotron::accumulate(
    string /*t3*/, float /*e3*/, 
    string /*t4*/, float e4)
{

  //if( (t1 == "e-" || t1 == "e+") && e1 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl
  //if( (t2 == "e-" || t2 == "e+") && e2 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl

  float f1 = 1.0f;
  float f2 = 3.0f;

  // 1/x^a suppression below threshold
  float thr = 1.0f;
  if( e4 < thr ) f2 *= pow(e4/thr, -2.0);

  return {f1,f2};
}
  


void Synchrotron::interact(
  string& /*t1*/, float& ux1, float& uy1, float& uz1,
  string& t2,     float& ux2, float& uy2, float& uz2) 
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
  float dx=0, dy, logchix, chi_x;
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
  // v2 with direct fitting formulas
  // not tested 

  //float a, n;
  //a  = 7.635/pow(chi_e, 0.996);
  //a *= 0.5 - 0.4995*tanh(0.1271*(chi_e - 21.93)); // a coeff
  //                                                  
  //n = 0.313*(1.0 + 3.2*exp(-pow(chi_e,0.4332))); // n coeff

  //float z = rand();
  //float chi_x2 = chi_e*pow(z/(1.0 + a*(1.0-z)), 1.0/n);

  ////std::cout << "chie:" << chi_e << " chi_x:" << chi_x << " chi_x2:" << chi_x2 << std::endl;
  //chi_x = chi_x2;

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
  // we assume that cooling follows du/dt = -P_0 u^2, then the solution is
  // u(t) = u_0/(1 - u_0*t/t_rad)
  // Here, t/t_rad = \Delta t/t_rad is suppliad via the wtar2wini parameter.

  //float A = 1.0f/(1.0f + z0*wtar2wini); // A = 1/(1 + u_0*(\Delta t/t_rad))
  //float facc = z0*(1.0f - A)/x; // f_acc = u_0*(1-A)/x_syn
  //facc = std::max(facc, 1.0f); // cap at 1

  // accumulated losses with a cap at 1e-3 change of an individual mom. component
  //ux1 -= std::min(facc*ux2, 0.999f*ux1);
  //uy1 -= std::min(facc*uy2, 0.999f*uy1);
  //uz1 -= std::min(facc*uz2, 0.999f*uz1);

  //wtar2wini = facc; // insert accumulation factor back to feed it into photon weight adjustment in pairing routines
  //std::cout << " syn emit energy:" << x << " pair gamma: " << z0 << " f_acc:" << wtar2wini << "\n";



  //--------------------------------------------------
  // v2 sub-step cooling via the energy evolution equation
    
  float N1Q = wtar2wini; // pairing routine gives us the onebody normalization coefficient via this parameter

  // radiative cooling time for the process per Delta t
  // t_rad/Delta t = 3/2 N gam^2/chi^2 \approx 3/2 N b^-2  for synchrotron
  // t_rad/Delta t = 3/2 N gam^2/chi^2 \approx 3/2 N b^-2 (rg/Rcurv)^-2 \gamma^-2 for virtual curvature
  //
  // NOTE: For virtual curvature radiation (in 1D) we feed the virtual pitch angle \sin\theta = gamma (rg/Rcurv) 
  //       into the mechanism via the B_y component that then gets included in the value of \chi_e. 
  // NOTE2: This expression ignores the quantum Gaunt factor and assumes b << 1
  float t_rad = 1.5*N1Q*(gam*gam/chi_e/chi_e);   // TODO

  // Solution of the curvature cooling energy equation is
  // A = g(Delta t)/g_0 = ( 1 + 3 * g^3 (Delta t/t_curv) )^-1/3
  //                    = ( 1 + 3 * g (Delta t/t_rad) )^-1/3 
  // since our definition of t_rad includes g^-2 factor already for virtual curvature process
  float A = 1.0f/cbrt(1.0f + 3.0f*gam/t_rad);

  // Solution of the synchrotron cooling energy equation is
  // A = g(Delta t)/g_0 = ( 1 + g (Delta t/t_rad) )^-1
  //float A = 1.0/(1.0 + gam/t_rad);

  // cap at a min of 1 and max of 10^-4 change
  A = std::max(1.0e-4f, std::min(1.0f, A)); 

  ////New electron u assuming only 1 photon emission from the regular MC process:
  float ux1b = ux1 - sign(ux1)*abs(ux2);
  float uy1b = uy1 - sign(uy1)*abs(uy2);
  float uz1b = uz1 - sign(uz1)*abs(uz2);

  ////Momentum reduction based on the cooling in one dt (but at least 0.001 of the original)
  //NOTE: we treat momentum as energy here assuming g >> 1
  float ux1c = A*ux1;
  float uy1c = A*uy1;
  float uz1c = A*uz1;

  //Making sure the new particle momentum is at least as small as based on 1 photon emission.
  ux1 = sign(ux1)*std::min(abs(ux1b), abs(ux1c));
  uy1 = sign(uy1)*std::min(abs(uy1b), abs(uy1c));
  uz1 = sign(uz1)*std::min(abs(uz1b), abs(uz1c));

  float gam_new = sqrt(1.0 + ux1*ux1 + uy1*uy1 + uz1*uz1);

  // growth factor; the this is how many reactions we did in reality during the unresolved dt
  float fgrowth = (gam - gam_new)/x;
  wtar2wini = fgrowth;


  //--------------------------------------------------
  // v1 with four momenta and curvature radiation 

  //--------------------------------------------------
  //if (ux1>0){
  //std::cout << " old vel:" << ux1 << "\n";}
  //double Ax = 1.0/cbrt(1.0 + sign(ux1)*3.0*C_SYNC*ux1*ux1*ux1);
  //double Ay = 1.0/cbrt(1.0 + sign(uy1)*3.0*C_SYNC*uy1*uy1*uy1);
  //double Az = 1.0/cbrt(1.0 + sign(uz1)*3.0*C_SYNC*uz1*uz1*uz1);

  //std::cout << "A:"   << Ax 
  //          << " vs.:" << A 
  //          << " rat:" << Ax/A 
  //          << " gam:" << gam
  //          << " 1/trad:" << 1.0/t_rad
  //          << " CSYN:" << C_SYNC*ux1*ux1
  //          << " 1+gam/trad:" << 1.0 + 3.0*gam/t_rad
  //          << "\n";


  //double ene_old = sqrt(ux1*ux1+uy1*uy1+uz1*uz1);

  ////New u assuming only 1 photon emission:
  //double ux1b = ux1 - sign(ux1)*std::min(double(abs(ux2)), abs(0.999*ux1));
  //double uy1b = uy1 - sign(uy1)*std::min(double(abs(uy2)), abs(0.999*uy1));
  //double uz1b = uz1 - sign(uz1)*std::min(double(abs(uz2)), abs(0.999*uz1));

  //////Momentum reduction based on the cooling in one dt (but at least 0.001 of the original)
  //double ux1c = std::max(Ax,0.001)*ux1;
  //double uy1c = std::max(Ay,0.001)*uy1;
  //double uz1c = std::max(Az,0.001)*uz1;

  ////Making sure the new particle momentum is at least as small as based on 1 photon emission.
  //ux1 = sign(ux1)*std::min(abs(ux1b),abs(ux1c));
  //uy1 = sign(uy1)*std::min(abs(uy1b),abs(uy1c));
  //uz1 = sign(uz1)*std::min(abs(uz1b),abs(uz1c));

  //double ene_new = sqrt(ux1*ux1 + uy1*uy1 + uz1*uz1);
  // //if (ux1>0){
  // //std::cout << " new vel:" << ux1 << "\n";}
  // float facc = (ene_old-ene_new)/x;
  // wtar2wini = facc;

  //wtar2wini = 1.0;

  //std::cout << " syn emit energy:" << x << " pair gamma: " << z0 << " f_acc:" << wtar2wini << "\n";

#ifdef DEBUG 

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

    std::cerr << "ERROR QNTM-SYNCH:" << std::endl;
    for(size_t i = 0; i < 5; i++) { std::cerr << i << " " << ts[i] << std::endl; }

    std::cerr << "logchix:" << logchix << " chi_x:" << chi_x <<  " chi_e:" <<  chi_e << std::endl;
    std::cerr << "zv" <<  zv << std::endl;
    std::cerr << "x" <<  x << " z0:" << z0 << std::endl;

    std::cerr << "ux1" <<  ux1 << " uy1:" << uy1 << " uz1:" << uz1 << std::endl;
    std::cerr << "ux2" <<  ux2 << " uy2:" << uy2 << " uz2:" << uz2 << std::endl;

    std::cerr << "i:" <<  i << " dx:" << dx << std::endl;
    std::cerr << "j:" <<  j << " dy:" << dy << std::endl;
  }


#endif

  return;
}


} // end of namespace qed

#pragma GCC diagnostic pop
