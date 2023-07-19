#pragma once

#include "interaction.h"


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
    cross_section = 0.35;  //0.68; //0.75; // 2*(3/8)*sigma_T
  }

  tuple<float_p, float_p> get_minmax_ene( string t1, string t2, double ene) override;

  float_p comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2) override;

  //tuple<
  //  string, float_p, float_p, float_p,
  //  string, float_p, float_p, float_p>
  //    interact(
  //      string t1, float_p ux1, float_p uy1, float_p uz1,
  //      string t2, float_p ux2, float_p uy2, float_p uz2) override;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override;


}; // end of PairAnn class


} // end of namespace qed
