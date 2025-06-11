#pragma once

#include "external/corgi/corgi.h"
#include "external/corgi/tile.h"
#include "tools/config_parser.h"

#include <string>
#include <unordered_map>

namespace emf2 {
/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : virtual public corgi::Tile<D> {

public:
  Tile(const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;
};

}  // namespace emf2
