#pragma once

#include "io/snapshots/fields.h"

namespace h5io {

/// IO object for storing (compressed) snapshots of basic emf2::YeeLattice quantities.
template<size_t D>
class FieldsWriter2 : public FieldsWriter<D> {

  using FieldsWriter<D>::FieldsWriter;

  /// read tile meshes into memory
  void read_tiles(corgi::Grid<D>& grid) override;
};

}  // end of namespace h5io
