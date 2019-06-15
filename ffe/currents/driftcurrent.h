#pragma once

#include "current.h"

namespace ffe {

/// Reduced system where only drift currents are solved
template<size_t D>
class DriftCurrent :
  public virtual Current<D>
{

  void interpolate_b(
    fields::YeeLattice& yee,
    ffe::ExtraLattice& tmp
    );

  void comp_drift_cur(Tile<D>& tile) override;

  void comp_parallel_cur(Tile<D>& tile) override;


};


} // end of namespace solve
