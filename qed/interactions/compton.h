#pragma once

#include "interaction.h"


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
  }

  tuple<float_p, float_p> get_minmax_ene( string t1, string t2 ) override;

  float_p comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2) override;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override;


}; // end of Compton class


} // end of namespace qed

