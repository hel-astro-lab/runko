#pragma once

#include "interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


class PairAnn :
  public Interaction
{
public:

  using Interaction::Interaction;

  tuple<double, double> get_minmax_ene( string t1, string t2 ) override;

  double comp_cross_section(
    string t1, double ux1, double uy1, double uz1,
    string t2, double ux2, double uy2, double uz2) override;

  //tuple<
  //  string, double, double, double,
  //  string, double, double, double>
  //    interact(
  //      string t1, double ux1, double uy1, double uz1,
  //      string t2, double ux2, double uy2, double uz2) override;

  void interact(
        string& t1, double& ux1, double& uy1, double& uz1,
        string& t2, double& ux2, double& uy2, double& uz2) override;


}; // end of PairAnn class


} // end of namespace qed
