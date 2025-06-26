#pragma once

#include "tyvi/mdgrid.h"
#include "tyvi/mdgrid_buffer.h"

#include <cstddef>
#include <experimental/mdspan>
#include <vector>

namespace emf2 {

/// Yee lattice of plasma quantities in tyvi::mdgrid continers.
struct [[nodiscard]] YeeLattice {
  using grid_extents_type = std::dextents<std::size_t, 3>;
  static constexpr auto vec_element =
    tyvi::mdgrid_element_descriptor<float> { .rank = 1, .dim = 3 };

  using vec_grid = tyvi::mdgrid<vec_element, grid_extents_type>;

  using host_vec_grid = tyvi::mdgrid_buffer<
    std::vector<vec_grid::value_type>,
    vec_grid::element_extents_type,
    vec_grid::element_layout_type,
    vec_grid::grid_extents_type,
    vec_grid::grid_layout_type>;

  /// Electric field
  vec_grid E;

  /// Magnetic field
  vec_grid B;

  /// Current
  vec_grid J;

  YeeLattice(std::size_t Nx, std::size_t Ny, std::size_t Nz);
};

}  // namespace emf2
