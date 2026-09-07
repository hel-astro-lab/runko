// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "entt/entt.hpp"
#include "runko/tools/config_parser.h"

#include <ranges>

namespace runko {
struct simulation_context {
  template<typename... Ts>
  std::ranges::view auto view_tiles();

  template<typename... Ts>
  std::ranges::view auto view_tiles() const;

  const toolbox::ConfigParser config;
  entt::registry tiles;
  using tile_id_type = entt::registry::entity_type;
};

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles()
{ return this->tiles.view<Ts...>().each(); }

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles() const
{ return this->tiles.view<const Ts...>().each(); }

}  // namespace runko
