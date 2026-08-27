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
  }

  // maximum of sigma_gg F/2 over all kinematics; 0.25563 at x_i x_j (1-mu) = 3.96
  const float cross_section = 0.256;

  tuple<float, float> get_minmax_ene( string t1, string t2, double ene) override final;

  pair_float comp_cross_section(
    string t1, float ux1, float uy1, float uz1,
    string t2, float ux2, float uy2, float uz2) override;

  void interact(
        string& t1, float& ux1, float& uy1, float& uz1,
        string& t2, float& ux2, float& uy2, float& uz2) override;


}; // end of PhotAnn class


} // end of namespace qed

