// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "tyvi/actions_ast.h"
#include "tyvi/actions_list.h"
#include "tyvi/execution.h"
#include "tyvi/sstd.h"

#include <format>
#include <ranges>
#include <source_location>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

namespace runko {

/// Parses arguments T... from given sexpr and returns a sender representing them.
///
/// Throws if a argument is missing or there is more arguments than sizeof...(T).
template<typename... T>
auto
  parse_atom_args(
    const tyvi::actions::sexpr& args,
    const std::source_location location = std::source_location::current())
{
  try {
    auto v       = tyvi::actions::list_view(args);
    const auto N = std::ranges::distance(v);
    if(N != sizeof...(T)) {
      throw std::runtime_error { std::format(
        "invalid number of arguments (expecting {}, got {})",
        sizeof...(T),
        N) };
    }

    auto atoms = std::vector<tyvi::actions::atom>();
    for(const auto& [n, s]: std::views::enumerate(v)) {
      if(not std::holds_alternative<tyvi::actions::atom>(s)) {
        throw std::runtime_error { std::format("argument {} is not a atom", n) };
      }

      atoms.push_back(std::get<tyvi::actions::atom>(s));
    }

    auto handle_atom = [&]<std::size_t n, typename U>() {
      if(auto x = tyvi::actions::atom_cast<U>(atoms[n])) {
        return std::move(x).value();
      }
      throw std::runtime_error { std::format(
        "atom at argument {} does not hold object of type: '{}' (note that the name is "
        "from typeid(T).name() and thus might hold compiler specific name for the "
        "type).",
        n,
        typeid(U).name()) };
    };

    return [&]<std::size_t... I>(std::index_sequence<I...>) {
      return tyvi::exec::just(handle_atom.template operator()<I, T>()...);
    }(std::make_index_sequence<sizeof...(T)>());

  } catch(const std::exception& e) {
    throw std::runtime_error {
      std::format("{}: argument parsing error: {}", location.function_name(), e.what())
    };
  }
}


}  // namespace runko
