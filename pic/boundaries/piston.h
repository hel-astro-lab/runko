#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"

namespace pic {

/// Reflecting piston wall
template<size_t D>
class Piston
{

  public:

  Piston() = default;

  double walloc = 0.0;     // x location of the wall
  double gammawall = 1.0;  // gamma of the moving wall
  double betawall  = 0.0;  // beta of the moving wall

  /// \brief interpolate electromagnetic fields to particle locations
  void solve(pic::Tile<D>&  /*tile*/);
    
  /// \brief apply conducting boundaries behind piston head
  void field_bc(pic::Tile<D>& );

  /// Small current deposition routine for individual particles 
  void zigzag(pic::Tile<D>& /*tile*/, 
      double x2glob, double y2glob, double z2glob, 
      double x1, double y1, double z1, 
      double q);
};


} // end of namespace pic
