#pragma once

#include "core/qed/interactions/interaction.h"


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
  const float cross_section = 0.256; // 1.37*(3/8)*sigma_T 

  tuple<float, float> get_minmax_ene( string t1, string t2, double ene) override final;

  pair_float comp_cross_section(
    string t1, float ux1, float uy1, float uz1,
    string t2, float ux2, float uy2, float uz2) override;

  void interact(
        string& t1, float& ux1, float& uy1, float& uz1,
        string& t2, float& ux2, float& uy2, float& uz2) override;


}; // end of PhotAnn class


} // end of namespace qed

