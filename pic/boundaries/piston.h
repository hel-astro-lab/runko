#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"

namespace pic {

/// Reflecting piston wall
template<size_t D>
class Piston
{

  public:

  Piston() {};

  double walloc;    // x location of the wall
  double gammawall; // gamma of the moving wall
  double betawall;  // beta of the moving wall

  /// \brief interpolate electromagnetic fields to particle locations
  void solve(pic::Tile<D>& );

  /// Small current deposition routine for individual particles 
  void zigzag(pic::Tile<D>&, 
      double x0, double y0, double z0, 
      double x1, double y1, double z1, 
      double q);
};


} // end of namespace pic
