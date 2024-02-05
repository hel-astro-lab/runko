#include <iostream>
#include <cmath>
#include <cassert>
#include "synchrotron.h"
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


tuple<float_p, float_p> Synchrotron::get_minmax_ene( string t1, double ene)
{

  // TODO
    
  return {0.0f, INF};
}


float_p Synchrotron::comp_optical_depth(
    string t1, 
    float_p ux1, float_p uy1, float_p uz1,
    float_p ex,  float_p ey,  float_p ez,
    float_p bx,  float_p by,  float_p bz)
{

  Vec3<float_p> zv;
  if( (t1 == "e-" || t1 == "e+") ) {
    zv.set(ux1, uy1, uz1); 
  } else { // incompatible particle types
    assert(false);
  }

  float_p gam = sqrt(1.0 + dot(zv,zv));   // prtcl gamma
  Vec3<float_p> beta = zv/gam;            // prtcl 3-velocity \beta

  // TODO
  float_p s0   = 1.0f;

  return s0;
}


//tuple<float_p, float_p> Synchrotron::accumulate(
//    string t1, float_p e1, string t2, float_p e2)
//{
//
//  //if( (t1 == "e-" || t1 == "e+") && e1 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl
//  //if( (t2 == "e-" || t2 == "e+") && e2 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl
//
//  // TODO
//
//  float_p f1 = 1.0f;
//  float_p f2 = 1.0f;
//
//  return {f1,f2};
//}
  


void Synchrotron::interact(
  string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
  string& t2, float_p& ux2, float_p& uy2, float_p& uz2) 
{

  //--------------------------------------------------
  Vec3<float_p> zv, xv;
  if( (t1 == "e-" || t1 == "e+") && (t2 == "ph") ) {
    zv.set(ux1, uy1, uz1); 
    xv.set(ux2, uy2, uz2); 
  } else if( (t2 == "e-" || t2 == "e+") && (t1 == "ph") ) {
    zv.set(ux2, uy2, uz2); 
    xv.set(ux1, uy1, uz1); 
  } else { // incompatible particle types
    assert(false);
  }

  //--------------------------------------------------

  return;
}


} // end of namespace qed
