// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "entt/entt.hpp"
#include "runko/tools/config_parser.h"
#include "runko/tools/timer.h"

#include <concepts>
#include <functional>
#include <optional>
#include <ranges>
#include <type_traits>
#include <typeinfo>

namespace runko {

template<typename F, typename T>
concept compatible_config_mapper =
  std::invocable<F, const toolbox::ConfigParser&> and
  std::convertible_to<std::invoke_result_t<F, const toolbox::ConfigParser&>, T>;

/// Simulation context is made of immutable config and entt::registry of tiles.
struct simulation_context {
  template<typename... Ts>
  std::ranges::view auto view_tiles();

  template<typename... Ts>
  std::ranges::view auto view_tiles() const;

  template<typename T>
  bool has_context_variable() const;

  /// Returns a reference to context variable T.
  ///
  /// Throws exception if the context variable does not exists.
  template<typename T>
  T& get_context_variable();

  /// Returns a reference to context variable T.
  ///
  /// Throws exception if the context variable does not exists.
  template<typename T>
  const T& get_context_variable() const;

  /// Set context variable using given value.
  ///
  /// If the context variable has already been set, does nothing.
  /// Returns true if the parameter is set and false if it was already set.
  template<typename T>
  bool set_context_variable(T&&);

  /// Set context variable using given mapper.
  ///
  /// If the context variable has already been set, does nothing.
  /// Returns true if the parameter is set and false if it was already set.
  template<typename T, compatible_config_mapper<T> F>
  bool set_context_variable(F&&);

  /// Returns a reference to context variable T.
  ///
  /// If context variable does not exists, potentially set it using given mapper.
  template<typename T, compatible_config_mapper<T> F>
  T& get_n_set_context_variable(F&&);

  const toolbox::ConfigParser config;
  entt::registry tiles;
  using tile_id_type = entt::registry::entity_type;

  runko::duration_container durations;
};

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles()
{ return this->tiles.view<Ts...>().each(); }

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles() const
{ return this->tiles.view<const Ts...>().each(); }

template<typename T>
bool
  simulation_context::has_context_variable() const
{ return this->tiles.ctx().contains<T>(); }

template<typename T>
T&
  simulation_context::get_context_variable()
{
  if(not this->tiles.ctx().contains<T>()) {
    throw std::runtime_error { std::format(
      "Config of a simulation_context does not contain: {}",
      typeid(T).name()) };
  }

  return this->tiles.ctx().get<T>();
}

template<typename T>
const T&
  simulation_context::get_context_variable() const
{
  if(not this->tiles.ctx().contains<T>()) {
    throw std::runtime_error { std::format(
      "Config of a simulation_context does not contain: {}",
      typeid(T).name()) };
  }

  return this->tiles.ctx().get<T>();
}

template<typename T>
bool
  simulation_context::set_context_variable(T&& x)
{
  using U = std::remove_cvref_t<T>;
  if(this->tiles.ctx().contains<U>()) { return false; }
  this->tiles.ctx().emplace<U>(std::forward<T>(x));
  return true;
}


template<typename T, compatible_config_mapper<T> F>
bool
  simulation_context::set_context_variable(F&& mapper)
{
  if(this->tiles.ctx().contains<T>()) { return false; }
  this->tiles.ctx().emplace<T>(std::invoke(std::forward<F>(mapper), this->config));
  return true;
}

template<typename T, compatible_config_mapper<T> F>
T&
  simulation_context::get_n_set_context_variable(F&& mapper)
{
  if(not this->tiles.ctx().contains<T>()) {
    this->tiles.ctx().emplace<T>(std::invoke(std::forward<F>(mapper), this->config));
  }

  return this->tiles.ctx().get<T>();
}
}  // namespace runko
