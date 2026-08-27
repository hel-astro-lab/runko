#pragma once

#include "core/qed/interactions/interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class PairAnn :
  public Interaction
{
public:

  PairAnn(string t1, string t2) : 
    Interaction(t1, t2)
  {
    name = "pair-ann";
  }

  // maximum of sigma_ann F/2 over all kinematics; 0.20798 at gamma_cm = 1.1354
  const float cross_section = 0.208;

  tuple<float, float> get_minmax_ene( string t1, string t2, double ene) override final;

  pair_float comp_cross_section(
    string t1, float ux1, float uy1, float uz1,
    string t2, float ux2, float uy2, float uz2) override;

  //tuple<
  //  string, float, float, float,
  //  string, float, float, float>
  //    interact(
  //      string t1, float ux1, float uy1, float uz1,
  //      string t2, float ux2, float uy2, float uz2) override;

  void interact(
        string& t1, float& ux1, float& uy1, float& uz1,
        string& t2, float& ux2, float& uy2, float& uz2) override;


}; // end of PairAnn class


} // end of namespace qed
