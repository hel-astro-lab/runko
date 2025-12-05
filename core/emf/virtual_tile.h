#pragma once

#include "core/emf/common.h"
#include "core/emf/yee_lattice.h"
#include "corgi/tile.h"
#include "mpi4cpp/mpi.h"
#include "tools/config_parser.h"
#include "tools/hollow_grid.h"

#include <array>
#include <cstddef>

namespace emf {

template<std::size_t D>
class Tile;

/*! \brief Virtual tile specialization of the emf tile.

 * Hollow version of the emf tile, i.e. stores only the halo regions.
 */
template<std::size_t D>
class VirtualTile : virtual public corgi::Tile<D> {
public:
  // FIXME: unify this to common source for emf::[Virtual]Tile.
  static_assert(D == 3);

  friend class Tile<D>;

  using value_type = emf::YeeLattice::value_type;

  using hollow_grid_EB = toolbox::hollow_grid<value_type, D, halo_size>;
  using hollow_grid_J  = toolbox::hollow_grid<value_type, D, 2 * halo_size>;

private:
  std::array<std::size_t, 3> extents_wout_halo_;

  hollow_grid_EB E_, B_;
  hollow_grid_J J_;

public:
  /// Construct VirtualTile by deducing extents from given tile grid index and config.
  ///
  /// Config has to contain values for:
  ///
  /// `N{x,y,z}`: tile grid extents
  /// `N{x,y,z}Mesh`: extents of mesh/grid in each tile
  /// `{x,y,z}min: minimum coordinate values
  /// `d{x,y,z}`: coordinate distance between mesh/grid points
  explicit VirtualTile(
    std::array<std::size_t, 3> tile_indices,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~VirtualTile() = default;

  /// Size of the non-halo yee lattice.
  std::array<std::size_t, 3> extents_wout_halo() const;

  /// Size of the yee lattice with halo.
  std::array<std::size_t, 3> extents_with_halo() const;

  /// Receive field specified with mode (see comm_mode).
  std::vector<mpi4cpp::mpi::request>
    recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag)
      override;
};

}  // namespace emf
