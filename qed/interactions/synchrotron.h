#pragma once

#include "interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class Synchrotron :
  public SingleInteraction
{

public:

  Synchrotron(string t1) : 
    SingleInteraction(t1)
  {
    name = "synchrotron";
  }

  pair_float get_minmax_ene( string t1, double ene) override final;

  pair_float comp_optical_depth( string t1, float_p ux1, float_p uy1, float_p uz1 ) override;

  //pair_float accumulate(string t1, float_p e1 ) override;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override;


}; // end of Synchrotron class


} // end of namespace qed
