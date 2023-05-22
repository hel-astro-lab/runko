#include <iostream>
#include <cmath>
#include <cassert>
#include "pair_ann.h"
#include "../../tools/vector.h"


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


tuple<float_p, float_p> Compton::get_minmax_ene( string t1, string t2)
{
  return {0.0, 3.0};
}

float_p Compton::comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2)
{

  return s0*fkin;
}

  
void Compton::interact(
  string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
  string& t2, float_p& ux2, float_p& uy2, float_p& uz2) 
{
  Vec3<float_p> zmvec(ux1, uy1, uz1);
  Vec3<float_p> zpvec(ux2, uy2, uz2);









  t1 = "ph";
  ux1 = xpp0(1);
  uy1 = xpp0(2);
  uz1 = xpp0(3);

  t2 = "ph";
  ux2 = xpp1(1);
  uy2 = xpp1(2);
  uz2 = xpp1(3);

  return;
}


} // end of namespace qed
