#include <cmath>
#include "pair_ann.h"


namespace qed {
  using std::string;
  using std::tuple;

tuple<double, double> PairAnn::get_minmax_ene( string t1, string t2)
{
  return {0.0, 3.0};
}


double PairAnn::comp_cross_section(
    string t1, double ux1, double uy1, double uz1,
    string t2, double ux2, double uy2, double uz2)
{
  return 1.0;
}


//tuple<
//    string, double, double, double,
//    string, double, double, double>
//PairAnn::interact(
//  string t1, double ux1, double uy1, double uz1,
//  string t2, double ux2, double uy2, double uz2) 
//{
//
//  auto ux3 = ux1;
//  auto uy3 = uy1;
//  auto uz3 = uz1;
//
//  auto ux4 = ux2;
//  auto uy4 = uy2;
//  auto uz4 = uz2;
//
//  auto t3 = t1;
//  auto t4 = t2;
//
//  //      particle 1          particle 2
//  return {t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4};
//}
  
void PairAnn::interact(
  string& t1, double& ux1, double& uy1, double& uz1,
  string& t2, double& ux2, double& uy2, double& uz2) 
{
  ux1 = 1.0;
  uy1 = 1.0;
  uz1 = 1.0;

  ux2 = 1.0;
  uy2 = 1.0;
  uz2 = 1.0;

  return;
}


} // end of namespace qed
