#pragma once

#include "core/emf2/tile.h"
#include "core/pic2/particle.h"
#include "external/corgi/corgi.h"
#include "external/corgi/tile.h"
#include "tools/config_parser.h"

#include <array>
#include <cstddef>
#include <string>

namespace pic2 {

namespace mpi = mpi4cpp::mpi;

/*! \brief PiC v2 tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf2::Tile
 */

template<std::size_t D>
class Tile : virtual public emf2::Tile<D>, virtual public corgi::Tile<D> {

  static_assert(D == 3);

  std::vector<ParticleContainer<D>> particle_buffs_;

  std::size_t particles_per_cell_;

public:
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_idx,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;
};


}  // namespace pic2
