#pragma once

#include "external/corgi/corgi.h"
#include "external/corgi/tile.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <cstddef>
#include <experimental/mdspan>
#include <functional>
#include <string>
#include <tuple>
#include <unordered_map>

namespace emf2 {

/// Yee lattice of plasma quantities
struct [[nodiscard]] YeeLattice {
  using GridExtents = std::dextents<std::size_t, 3>;
  static constexpr auto VecElement =
    tyvi::mdgrid_element_descriptor<float> { .rank = 1, .dim = 3 };
  static constexpr auto ScalarElement =
    tyvi::mdgrid_element_descriptor<float> { .rank = 0, .dim = 3 };

  using VecGrid    = tyvi::mdgrid<VecElement, GridExtents>;
  using ScalarGrid = tyvi::mdgrid<ScalarElement, GridExtents>;

  /// Electric fiel
  VecGrid E;

  /// Magnetic field
  VecGrid B;

  /// Charge density
  ScalarGrid rho;

  /// Current vector
  VecGrid J;

  YeeLattice(std::size_t Nx, std::size_t Ny, std::size_t Nz);
};

/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in
 * staggered Yee lattice.
 */
template<std::size_t D>
class Tile : virtual public corgi::Tile<D> {

  static_assert(D <= 3);

protected:
  std::array<std::size_t, 3> yee_lattice_extents_;

private:
  YeeLattice yee_lattice_;

  struct with_halo_tag {};
  struct wout_halo_tag {};
  static constexpr with_halo_tag with_halo {};
  static constexpr wout_halo_tag wout_halo {};

  /// Returns tuple of mdspans to yee_lattice_ staging buffer {E, B, rho, J}.
  [[nodiscard]] auto yee_lattice_staging_mds(with_halo_tag)
  {
    return std::tuple { yee_lattice_.E.staging_mds(),
                        yee_lattice_.B.staging_mds(),
                        yee_lattice_.rho.staging_mds(),
                        yee_lattice_.J.staging_mds() };
  }

  [[nodiscard]] auto yee_lattice_staging_mds(wout_halo_tag)
  {
    auto [Emds, Bmds, rhomds, Jmds] = yee_lattice_staging_mds(with_halo);

    const auto x = std::tuple { halo_length, halo_length + yee_lattice_extents_[0] };
    const auto y = std::tuple { halo_length, halo_length + yee_lattice_extents_[1] };
    const auto z = std::tuple { halo_length, halo_length + yee_lattice_extents_[2] };

    return std::tuple { std::submdspan(std::move(Emds), x, y, z),
                        std::submdspan(std::move(Bmds), x, y, z),
                        std::submdspan(std::move(rhomds), x, y, z),
                        std::submdspan(std::move(Jmds), x, y, z) };
  }

  double cfl_;

public:
  static constexpr auto halo_length = 3;

  Tile(const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;

  using scalar_field = std::function<double(double, double, double)>;
  using vector_field = std::function<std::array<double, 3>(double, double, double)>;

  /// Precondition: corgi::Tile<D>::{mins,maxs} are set.
  void set_fields(vector_field E, vector_field B, scalar_field rho, vector_field J);
};

}  // namespace emf2
