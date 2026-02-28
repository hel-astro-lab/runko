#pragma once


#include "core/mdgrid_common.h"
#include "thrust/device_vector.h"
#include "tyvi/mdgrid.h"


namespace vlv {


// input arguments used fed to the class constructor
struct VlasovGridArgs {
  std::size_t N;
  double charge, mass;
};



class [[nodiscard]] VlasovGrid {
public:
  /// The type in which pos and vel are stored in.
  using value_type = float;

public:
  explicit VlasovGrid(VlasovGridArgs);



}





}  // namespace vlv
