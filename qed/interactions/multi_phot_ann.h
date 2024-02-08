#pragma once

#include "interaction.h"


namespace qed {
  using std::string;
  using std::tuple;


 // Multi-photon Breit-Wheeler process 
 // also known as magnetic single-photon pair creation
 //
 // Refs Duclous et al. 2011 implementation
 //
class MultiPhotAnn :
  public Interaction
{

private:

  //float CHIE[32] = 

  //float XI[32][33] = {


public:

  MultiPhotAnn(string t1) : 
    Interaction(t1, "")
  {
    name = "multi-phot-ann";
    interaction_order = 1; // set as single prtcl interaction
  }

  double B_QED = 1.0; // critical QED field strength in code units

  // internal storage for the dimensionless photon quantum parameter
  // updated in comp_chi()
  float_p chi_x = 0.0; 

  // NOTE no override since input arguments are different
  pair_float get_minmax_ene( string t1, string t2, double ene) override final;

  // calculate quantum parameter
  float_p comp_chi( 
      float_p ux1, float_p uy1, float_p uz1,
      float_p ex,  float_p ey,  float_p ez,
      float_p bx,  float_p by,  float_p bz);

  // calculate optical depth for the process 
  float_p comp_optical_depth( string t1, 
      float_p ux1, float_p uy1, float_p uz1,
      float_p ex,  float_p ey,  float_p ez,
      float_p bx,  float_p by,  float_p bz
      ) override final;

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override final;

  //pair_float accumulate(string t1, float_p e1 ) override;

}; // end of MultiPhotAnn class

} // end of namespace qed

