#pragma once

#include "core/emf2/yee_lattice.h"
#include "corgi/corgi.h"
#include "corgi/tile.h"
#include "tools/config_parser.h"

#include <array>
#include <cstddef>
#include <experimental/mdspan>
#include <functional>
#include <tuple>

namespace emf2 {

enum class FieldPropagator { FDTD2 };

/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : virtual public corgi::Tile<D> {

  static_assert(D == 3);

protected:
  YeeLattice yee_lattice_;
  double cfl_;
  FieldPropagator field_propagator_;

public:
  static constexpr auto halo_size = 3;

  /// Construct Tile by deducing extents from given tile grid index and config.
  ///
  /// Config has to contain values for:
  ///
  /// `cfl`: (FIXME: add description)
  /// `N{x,y,z}`: tile grid extents
  /// `N{x,y,z}Mesh`: extents of mesh/grid in each tile
  /// `{x,y,z}min: minimum coordinate values
  /// `d{x,y,z}`: coordinate distance between mesh/grid points
  /// `field_propagator`: scheme to propagate E and B fields.
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_indices,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;

  using vector_field_function =
    std::function<std::tuple<double, double, double>(double, double, double)>;

  /// Set E, B and J fields in non-halo regions.
  void
    set_EBJ(vector_field_function E, vector_field_function B, vector_field_function J);

  /// Get E, B and J fields in non-halo regions.
  YeeLattice::YeeLatticeHostCopy get_EBJ();

  /// Size of the non-halo yee lattice.
  std::array<std::size_t, 3> extents_wout_halo() const;

  /// Advance B by half time step using scheme from configuration in non-halo region.
  void push_half_b();

  /// Advance E by full time step using scheme from configuration in non-halo region.
  ///
  /// Does not add contributions from the current.
  void push_e();

  /// E += J
  void deposit_current();
};

}  // namespace emf2
