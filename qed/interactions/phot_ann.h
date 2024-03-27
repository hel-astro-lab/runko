#pragma once

#include "qed/interactions/interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class PhotAnn :
  public Interaction
{
public:

  PhotAnn(string t1, string t2) : 
    Interaction(t1, t2)
  {
    name = "phot-ann";
    //cross_section = 0.256; // 0.682; //0.51375; // 1.37*(3/8)*sigma_T // FIXME
    //                       // 0.25564 measured
  }

  // maximum cross section
  const float_p cross_section = 0.256; // 1.37*(3/8)*sigma_T 

  tuple<float_p, float_p> get_minmax_ene( string t1, string t2, double ene) override final;

  pair_float comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2) override;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override;


}; // end of PhotAnn class


} // end of namespace qed

