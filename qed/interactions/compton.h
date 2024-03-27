#pragma once

#include "qed/interactions/interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class Compton :
  public Interaction
{
public:

  Compton(string t1, string t2) : 
    Interaction(t1, t2)
  {
    name = "compton";
    //cross_section = 1.0;  // maximum cross section 
                          // for head-on collisions its x2 (in units of sigma_T)
  }

  //maximum cross section 
  const float cross_section = 1.0;

  bool no_electron_update = false; // do not update electrons
  bool no_photon_update   = false; // do not update photons

  double ming = 1.1;      // minimumjj electron energy to classify it "non-relativistic"
  double minx2z = 1.0e-2; // minimum ph energy needs to be > minx2z*gam 

  tuple<float, float> get_minmax_ene( string t1, string t2, double ene) override final;

  pair_float comp_cross_section(
    string t1, float ux1, float uy1, float uz1,
    string t2, float ux2, float uy2, float uz2) override;

  pair_float accumulate(string t1, float e1, string t2, float e2) override;

  void interact(
        string& t1, float& ux1, float& uy1, float& uz1,
        string& t2, float& ux2, float& uy2, float& uz2) override;


}; // end of Compton class


} // end of namespace qed

