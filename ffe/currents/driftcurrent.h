#pragma once

#include "current.h"

namespace ffe {

/// Reduced system where only drift currents are solved
template<size_t D>
class DriftCurrent :
  public Current<D>
{

  // base-class initialization needed to initialize bf meshes
  using Current<D>::Current;

  void interpolate_b(fields::YeeLattice& yee);

  void comp_drift_cur(Tile<D>& tile) override;

  void comp_parallel_cur(Tile<D>& tile) override;

  void limiter(Tile<D>& tile) override;


};


} // end of namespace solve
