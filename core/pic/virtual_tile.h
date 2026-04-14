#pragma once

#include "core/emf/virtual_tile.h"
#include "core/particles_common.h"
#include "core/pic/particle.h"
#include "corgi/tile.h"
#include "mpi4cpp/mpi.h"
#include "thrust/device_vector.h"
#include "tools/config_parser.h"

#include <array>
#include <vector>
#include <cstddef>

namespace pic {

template<std::size_t D>
class Tile;

/*! \brief Virtual tile specialization of the pic tile.
 *
 * Hollow version of the pic tile, i.e. stores only the particles in halo regions.
 */
template<std::size_t D>
class VirtualTile : virtual public emf::VirtualTile<D> {
public:
  static_assert(D == 3);

  friend class Tile<D>;

  using value_type = ParticleContainer::value_type;

private:
  std::size_t amount_of_ptypes_;
  std::vector<std::size_t> subregion_ends_;
  thrust::device_vector<runko::ParticleState<value_type>> buffer_;

public:
  /// Construct VirtualTile by deducing extents from given tile grid index and config.
  ///
  /// Config has to contain values for:
  ///
  /// `N{x,y,z}`: tile grid extents
  /// `N{x,y,z}Mesh`: extents of mesh/grid in each tile
  /// `{x,y,z}min: minimum coordinate values
  /// `d{x,y,z}`: coordinate distance between mesh/grid points
  ///
  /// Particle charges and masses as qx and mx where x is 0, 1, ...
  /// Note that particle charges and masses qx and mx are read in order: 0, 1, ...
  /// If a i:th mass and charge are missing, the search is stopped.
  explicit VirtualTile(
    std::array<std::size_t, 3> tile_indices,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~VirtualTile() = default;

  /// Receive field specified with mode (see comm_mode).
  std::vector<mpi4cpp::mpi::request>
    recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag)
      override;
};

}  // namespace pic
