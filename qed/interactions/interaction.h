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
  std::uniform_real_distribution<double> uni_dis;

public:
    
  double cross_section = 1.0; // maximum cross section (in units of sigma_T)

  string t1; // incident particle type
  string t2; // target particle type

  // constructor with incident/target types
  Interaction(string t1, string t2) :
    t1(t1),
    t2(t2),
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)       
  { }

  // minimum and maximum particle energies required to participate in the interaction
  virtual tuple<double, double> get_minmax_ene( string t1, string t2 ) { return {0.0, 1.0}; };

  // interaction cross section given incident/target particle four-velocities
  virtual double comp_cross_section(
    string t1, double ux1, double uy1, double uz1,
    string t2, double ux2, double uy2, double uz2)
    { return cross_section; }


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
  //  string, double, double, double,
  //  string, double, double, double>
  //    interact(
  //      string t1, double ux1, double uy1, double uz1,
  //      string t2, double ux2, double uy2, double uz2)
  //    { return 
  //      // particle 1        particle 2
  //      {t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2};
  //    }

  virtual void interact(
        string& t1, double& ux1, double& uy1, double& uz1,
        string& t2, double& ux2, double& uy2, double& uz2)
      { return; }

  // random numbers between [0, 1]
  double rand() { return uni_dis(gen); };


}; // end of class Interaction


} // end of namespace qed
