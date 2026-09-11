// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "entt/entt.hpp"
#include "runko/tools/config_parser.h"

#include <concepts>
#include <functional>
#include <optional>
#include <ranges>
#include <type_traits>

namespace runko {

template<typename F, typename T>
concept compatible_config_mapper =
  std::invocable<F, const toolbox::ConfigParser&> and
  std::convertible_to<std::invoke_result_t<F, const toolbox::ConfigParser&>, T>;

struct simulation_context {
  template<typename... Ts>
  std::ranges::view auto view_tiles();

  template<typename... Ts>
  std::ranges::view auto view_tiles() const;


  /// Returns a reference to context variable T.
  ///
  /// If context variable does not exists, potentially set it using given mapper.
  template<typename T, compatible_config_mapper<T> F>
  const T& get_n_set_config(F&&);

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


template<typename T, compatible_config_mapper<T> F>
const T&
  simulation_context::get_n_set_config(F&& mapper)
{
  if(not this->tiles.ctx().contains<T>()) {
    this->tiles.ctx().emplace<T>(std::invoke(std::forward<F>(mapper), this->config));
  }

  return this->tiles.ctx().get<T>();
}
}  // namespace runko
