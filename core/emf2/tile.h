#pragma once

#include "external/corgi/corgi.h"
#include "external/corgi/tile.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <cstddef>
#include <experimental/mdspan>
#include <string>
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
  using ScalarGrid = tyvi::mdgrid<VecElement, GridExtents>;

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

  double cfl_;

public:
  Tile(const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;
};

}  // namespace emf2
