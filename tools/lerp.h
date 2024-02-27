#pragma once

namespace toolbox {

// linear cubic interpolation routine
inline double lerp(
      double c000,
      double c100,
      double c010,
      double c110,
      double c001,
      double c101,
      double c011,
      double c111,
      double dx, double dy, double dz
      ) 
{
      double c00 = c000 * (1.0-dx) + c100 * dx;
      double c10 = c010 * (1.0-dx) + c110 * dx;
      double c0  = c00  * (1.0-dy) + c10  * dy;
      double c01 = c001 * (1.0-dx) + c101 * dx;
      double c11 = c011 * (1.0-dx) + c111 * dx;
      double c1  = c01  * (1.0-dy) + c11  * dy;
      double c   = c0   * (1.0-dz) + c1   * dz;
      return c;
}


}
