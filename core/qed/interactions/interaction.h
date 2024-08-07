#pragma once

#include <string>
#include <tuple>
#include <random>

#include "definitions.h"

// TODO turning compiler warnings off temporarily in this file since 
//      for symmetry, there are lots of unused variables in the qed API
// NOTE remember to remove pragma pop at the end of the file when done with these.
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wunused-but-set-variable"


namespace qed {
  using std::string;
  using std::tuple;



// Base class for a generic two-body QED interaction
//class SingleInteraction
//{
//private:
//
//  std::random_device rd;
//  std::mt19937 gen;
//  std::uniform_real_distribution<float> uni_dis;
//
//public:
//
//  using pair_float = std::tuple<float, float>;
//
//  string name = "GenSingleInt"; // interaction name
//
//  string t1; // incident particle type
//
//  // error tolerance
//  const float tol = 3.0e-4;
//
//  // constructor with incident/target types
//  SingleInteraction(string t1) :
//    t1(t1),
//    gen(42), // gen(rd() ) 
//    uni_dis(0.0, 1.0)       
//  { }
//
//  // minimum and maximum particle energies required to participate in the interaction
//  virtual pair_float get_minmax_ene( string t1, double ene) { return {0.0f, 1.0f}; };
//
//  // interaction cross section given incident particle's four-velocities
//  virtual pair_float comp_optical_depth(
//    string t1, float ux1, float uy1, float uz1)
//    { return {1.0f, 1.0f}; }
//
//  // interaction accumulation factor
//  //virtual pair_float accumulate(string t1, float e1) {return {1.0f,1.0f}; };
//
//  virtual void interact(
//        string& t1, float& ux1, float& uy1, float& uz1,
//        string& t2, float& ux2, float& uy2, float& uz2)
//      { return; }
//
//  // random numbers between [0, 1[
//  float rand() { return uni_dis(gen); };
//
//  // random numbers between [a, b[
//  float rand_ab(float a, float b) { 
//    float r = rand();
//    return a + (b-a)*r;
//  };
//
//
//}; // end of class SingleInteraction



// Base class for a generic two-body QED interaction
class Interaction
{
private:

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float> uni_dis;

public:

  // constants used in calculations
  // NOTE: both constants are normalized to unity because of how the external (python-defined)
  //       normalization of QED reactions is done. It assumes that the numerical reaction rates 
  //       are of form alpha/lamC d\tau/dt (and so the alpha/lamC is dealt there).
  const float alphaf = 1.0; //1.0/137.0; // fine structure constant
  const float lamC = 1.0f; // normalized reduced Compton wavelength (real value is considered in pic.py)


  int interaction_order = 2; // default interaction is 2-body

  using pair_float = std::tuple<float, float>;
    
  const float cross_section = 1.0; // maximum cross section (in units of sigma_T)

  string name = "GenTwoBodyInt"; // interaction name
             //
  string t1; // incident particle type
  string t2; // target particle type

  // use accumulation technique (calls accumulate() function to get facc
  bool do_accumulate = false;

  // enhance energy transfer by this factor 
  float wtar2wini = 1.0f;

  // error tolerance
  const float tol = 3.0e-4;

  // maximum Monte Carlo iterations of the differential cross section
  const int n_max = 10000;

  // constructor with incident/target types
  Interaction(string t1, string t2) :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0),
    t1(t1),
    t2(t2) 
  { }

  // minimum and maximum particle energies required to participate in the interaction
  virtual pair_float get_minmax_ene( string /*t1*/, string /*t2*/, double /*ene*/) { return {0.0f, 1.0f}; };

  // interaction cross section given incident/target particle four-velocities
  // NOTE: used in binary interactions
  virtual pair_float comp_cross_section(
    string /*t1*/, float /*ux1*/, float /*uy1*/, float /*uz1*/,
    string /*t2*/, float /*ux2*/, float /*uy2*/, float /*uz2*/)
    { return {cross_section, 1.0f}; }

  // NOTE: used in single interactions
  virtual float comp_optical_depth( 
      string /*t1*/, 
      float /*ux1*/, float /*uy1*/, float /*uz1*/,
      float /*ex*/,  float /*ey*/,  float /*ez*/,
      float /*bx*/,  float /*by*/,  float /*bz*/)
  { return 1.0f; };

  // interaction accumulation factor
  virtual pair_float accumulate(string /*t1*/, float /*e1*/, string /*t2*/, float /*e2*/) {return {1.0f,1.0f}; };

  // main interaction routine; 
  //
  // takes as input: 
  //  particle1: type, ux, uy, uz
  //  particle2: type, ux, uy, uz
  //
  //  and returns:
  //  particle3: type, ux, uy, uz
  //  particle4: type, ux, uy, uz
  //
  virtual void interact(
        string& /*t1*/, float& /*ux1*/, float& /*uy1*/, float& /*uz1*/,
        string& /*t2*/, float& /*ux2*/, float& /*uy2*/, float& /*uz2*/)
      { return; }

  // random numbers between [0, 1[
  float rand() { return uni_dis(gen); };

  // random numbers between [a, b[
  float rand_ab(float a, float b) { 
    float r = rand();
    return a + (b-a)*r;
  };

}; // end of class Interaction


} // end of namespace qed
    
#pragma GCC diagnostic pop
