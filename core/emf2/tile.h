#pragma once

#include "external/corgi/corgi.h"
#include "external/corgi/tile.h"
#include "tools/config_parser.h"
#include "core/emf2/yee_lattice.h"

#include <array>
#include <cstddef>
#include <experimental/mdspan>
#include <functional>
#include <tuple>

namespace emf2 {

/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : virtual public corgi::Tile<D> {

  static_assert(D == 3);

protected:
  std::array<std::size_t, 3> yee_lattice_extents_wout_halo_;

private:
  YeeLattice yee_lattice_;

  /// Returns tuple of mdspans to yee_lattice_ staging buffer {E, B, J}.
  [[nodiscard]] auto yee_lattice_staging_mds_with_halo()
  {
    return std::tuple { yee_lattice_.E.staging_mds(),
                        yee_lattice_.B.staging_mds(),
                        yee_lattice_.J.staging_mds() };
  }

  [[nodiscard]] auto yee_lattice_staging_mds_wout_halo()
  {
    auto [Emds, Bmds, Jmds] = yee_lattice_staging_mds_with_halo();

    const auto x =
      std::tuple { halo_length, halo_length + yee_lattice_extents_wout_halo_[0] };
    const auto y =
      std::tuple { halo_length, halo_length + yee_lattice_extents_wout_halo_[1] };
    const auto z =
      std::tuple { halo_length, halo_length + yee_lattice_extents_wout_halo_[2] };

    return std::tuple { std::submdspan(std::move(Emds), x, y, z),
                        std::submdspan(std::move(Bmds), x, y, z),
                        std::submdspan(std::move(Jmds), x, y, z) };
  }

  double cfl_;

public:
  static constexpr auto halo_length = 3;

  /// Construct Tile by deducing extents from given tile grid index and config.
  ///
  /// Config has to contain values for:
  ///
  /// `cfl`: (FIXME: add description)
  /// `N{x,y,z}`: tile grid extents
  /// `N{x,y,z}Mesh`: extents of mesh/grid in each tile
  /// `{x,y,z}min: minimum coordinate values
  /// `d{x,y,z}`: coordinate distance between mesh/grid points
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

  struct host_EBJ_grids {
    YeeLattice::host_vec_grid E;
    YeeLattice::host_vec_grid B;
    YeeLattice::host_vec_grid J;
  };

  /// Get E, B and J fields in non-halo regions.
  host_EBJ_grids get_EBJ();

  /// Size of the non-halo yee lattice.
  std::array<std::size_t, 3> extents_wout_halo() const;
};

}  // namespace emf2
