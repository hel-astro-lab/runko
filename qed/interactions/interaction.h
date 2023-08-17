#pragma once

#include <string>
#include <tuple>
#include <random>
#include "../../definitions.h"


namespace qed {
  using std::string;
  using std::tuple;


// Base class for two-body a generic QED interaction
class Interaction
{
private:

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float_p> uni_dis;

public:
    
  float_p cross_section = 1.0; // maximum cross section (in units of sigma_T)

  string name = "GenInteraction"; // interaction name
             //
  string t1; // incident particle type
  string t2; // target particle type

  // use accumulation technique (calls accumulate() function to get facc
  bool do_accumulate = false;

  // constructor with incident/target types
  Interaction(string t1, string t2) :
    t1(t1),
    t2(t2),
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)       
  { }

  // minimum and maximum particle energies required to participate in the interaction
  virtual tuple<float_p, float_p> get_minmax_ene( string t1, string t2, double ene) { return {0.0f, 1.0f}; };

  // interaction cross section given incident/target particle four-velocities
  using pair_float = std::tuple<float_p, float_p>;
  virtual pair_float comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2)
    { return {cross_section, 1.0}; }

  // interaction accumulation factor
  virtual pair_float accumulate(string t1, float_p e1, string t2, float_p e2) {return {1.0f,1.0f}; };

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
  //virtual tuple<
  //  string, float_p, float_p, float_p,
  //  string, float_p, float_p, float_p>
  //    interact(
  //      string t1, float_p ux1, float_p uy1, float_p uz1,
  //      string t2, float_p ux2, float_p uy2, float_p uz2)
  //    { return 
  //      // particle 1        particle 2
  //      {t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2};
  //    }

  virtual void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2)
      { return; }

  // random numbers between [0, 1[
  float_p rand() { return uni_dis(gen); };

  // random numbers between [a, b[
  float_p rand_ab(float_p a, float_p b) { 
    float_p r = rand();
    return a + (b-a)*r;
  };

}; // end of class Interaction


} // end of namespace qed
