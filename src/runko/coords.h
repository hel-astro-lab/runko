// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "runko/tools/vector.h"

#include <array>
#include <cstddef>

namespace runko {

/// Maps tile local coordinates of tile_idx to global coordinates.
///
/// If global coordinates are from [0, n_tiles[i] * N[i]) for each dimension i,
/// then tile local coordinates of a tile at tile_idx are the global coordinates
/// but shifted by n_cells[i] * tile_idx[i].
template<std::size_t rank>
struct [[nodiscard]] global_coordinates_closure {
  using vec = toolbox::VecD<double, rank>;
  vec n_cells;
  vec tile_idx;

  template<typename... I>
  [[nodiscard]]
  constexpr std::array<double, rank> operator()(I... idx) const
  {
    const auto coeff = [&]<std::size_t... J>(std::index_sequence<J...>) {
      return vec((static_cast<double>(idx) / static_cast<double>(n_cells[J]))...);
    }(std::make_index_sequence<rank>());
    return ((coeff + tile_idx) * n_cells).data;
  }

  /// Shorthand for operator()(0...).
  [[nodiscard]]
  constexpr std::array<double, rank> mins() const
  {
    return [this]<std::size_t... J>(std::index_sequence<J...>) {
      return (*this)((0 * J)...);
    }(std::make_index_sequence<rank>());
  }

  /// Shorthand for operator()(n_cells[i]...).
  [[nodiscard]]
  constexpr std::array<double, rank> maxs() const
  {
    return [this]<std::size_t... J>(std::index_sequence<J...>) {
      return (*this)((this->n_cells[J])...);
    }(std::make_index_sequence<rank>());
  }
};

/// Constructs global_coordinates_closure from n_cells_per_tiles and given index.
global_coordinates_closure<3>
  global_coordinates(const simulation_context&, const std::array<double, 3>&);

/// Shorthand for { (0...), (n_tiles[i] * n_cells_per_tile[i]... ) }.
template<std::size_t rank>
std::array<std::array<double, rank>, 2>
  global_coordinate_extents(const simulation_context& sim)
{

  const auto cells = toolbox::get_extent_list(sim.config, "n_cells_per_tile", rank);
  const auto tiles = toolbox::get_extent_list(sim.config, "n_tiles", rank);

  using arr = std::array<double, rank>;

  return [&]<std::size_t... J>(std::index_sequence<J...>) {
    return std::array { arr { static_cast<double>(0 * J)... },
                        arr { static_cast<double>(tiles[J] * cells[J])... } };
  }(std::make_index_sequence<rank>());
}
}  // namespace runko
