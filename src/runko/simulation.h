// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include "entt/entt.hpp"
#include "runko/tools/vector.h"

#include <array>
#include <cstddef>
#include <print>
#include <ranges>

namespace runko {

template<std::size_t rank>
using cartesian_index = toolbox::VecD<std::size_t, rank>;

struct simulation_context {
  void foobar() { std::println("foobar"); }

  template<typename... Ts>
  std::ranges::view auto view_tiles();

  entt::registry tiles;
};

template<std::size_t rank>
void add_cartesian_tiles(simulation_context&, std::array<std::size_t, rank> extents);

// Implementation:

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles()
{ return this->tiles.view<Ts...>().each(); }

template<std::size_t rank>
void
  add_cartesian_tiles(
    simulation_context& sim,
    [[maybe_unused]] const std::array<std::size_t, rank> extents)
{

  using index_type = runko::cartesian_index<3>;

  for(const auto idx: tyvi::sstd::index_space_view(
        std::layout_right::mapping(std::dextents<std::size_t, 3>(3, 2, 4)))) {
    sim.tiles.emplace<index_type>(sim.tiles.create(), index_type { idx });
  }
}

}  // namespace runko
