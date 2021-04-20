#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"

namespace pic {

/// Reflecting piston wall
template<size_t D>
class PistonZdir
{

  public:

  PistonZdir() = default;

  double wallocz = 0.0; // z location of the wall
  int wdir = 1; // direction towards positive z

  /// \brief interpolate electromagnetic fields to particle locations
  void solve(pic::Tile<D>& tile);
    
  /// \brief apply conducting boundaries behind piston head
  void field_bc(pic::Tile<D>& );

  /// Small current deposition routine for individual particles 
  void zigzag(pic::Tile<D>& tile, 
      double x1, double y1, double z1, 
      double x2, double y2, double z2, 
      double q);
};


} // end of namespace pic
