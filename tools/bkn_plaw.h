#pragma once

#include <cmath>


namespace toolbox {

inline float bkn_plaw(
    float x,    // value
    float A,    // norm
    float x_b,  // x-location of break
    float a_1,  // index of first slope
    float a_2,  // index of second slope
    float delta)// smoothing lenght
{
  a_1 *= -1.0f;
  a_2 *= -1.0f;
  return A*( pow(x/x_b, -a_1)*pow(0.5f*(1.0f + pow(x/x_b, 1.0f/delta)), (a_1-a_2)*delta));
}


} // end of ns
