// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/vector.h"

#include <array>

namespace runko {

/// Returns a lambda that maps tile coordinates of give indices to global coordinates.
///
/// For tile index (i, j, k), tile coordinates are defined to be (0, 0, 0)
/// at the cell indices (0, 0, 0) of the tile.
///
/// For example, global_coordinate_map(..., {1, 2, 3})(0, 0, 0) == {Dx, 2 * Dy, 3 * Dz}.
/// where Dx, Dy and Dz are the dimensions of a one tile.
///
/// The lambda will return the coordinates as std::array<double, 3>
/// and it can be invoked even with fractional indices,
/// i.e. global_coordinates(...)(0.5, 0.5, 0.5) maps to middle of the cell at (0, 0, 0).
auto global_coordinates(const simulation_context&, const std::array<double, 3>&);

// Implementation:

auto
  global_coordinates(const simulation_context& sim, const std::array<double, 3>& idx)
{
  const auto cells  = toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);
  const auto vcells = toolbox::Vec3<double>(cells[0], cells[1], cells[2]);
  const auto vidx   = toolbox::Vec3<double>(idx);

  return [vcells, vidx](const auto i, const auto j, const auto k) {
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(vcells[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(vcells[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(vcells[2]);
    const auto coeff   = toolbox::Vec3<double>(x_coeff, y_coeff, z_coeff);
    return ((coeff + vidx) * vcells).data;
  };
}

}  // namespace runko
