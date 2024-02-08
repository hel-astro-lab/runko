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
    
  // TODO add numerical prefactor

  float_p chix = comp_chi( ux1,uy1,uz1, ex,ey,ez, bx,by,bz); // dimensionless quantum parameter \chi_x

  // TODO

  return 1.0;
}



void MultiPhotAnn::interact(
  string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
  string& t2, float_p& ux2, float_p& uy2, float_p& uz2) 
{

  //--------------------------------------------------
  Vec3<float_p> xv, zm, zp;
  xv.set(ux1, uy1, uz1);  // 4-velocity 3-vector
  float_p x0 = norm(xv); // prtcl gamma
                                          
  //--------------------------------------------------
  // draw electron chi when photon chi is known
  
  // TODO

  return;
}


} // end of namespace qed
