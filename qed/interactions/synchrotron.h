#pragma once

#include "interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class Synchrotron :
  public Interaction
{

public:

  const int interaction_order = 1; // single prtcl interaction

  Synchrotron(string t1) : 
    Interaction(t1, "")
  {
    name = "synchrotron";
  }

  // NOTE no override since input arguments are different
  pair_float get_minmax_ene( string t1, double ene);

  float_p comp_optical_depth( string t1, 
      float_p ux1, float_p uy1, float_p uz1,
      float_p ex,  float_p ey,  float_p ez,
      float_p bx,  float_p by,  float_p bz
      ) override final;

  //pair_float accumulate(string t1, float_p e1 ) override;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override final;


}; // end of Synchrotron class


} // end of namespace qed

