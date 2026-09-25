// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pybind11/pybind11.h"

#include <concepts>
#include <format>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace toolbox {

class [[nodiscard]] ConfigParser {
  struct none_tag_type {};
  using value_type = std::variant<
    none_tag_type,
    bool,
    std::string,
    std::ptrdiff_t,
    double,
    std::vector<double>,
    std::vector<std::ptrdiff_t>,
    std::vector<std::string>>;
  std::unordered_map<std::string, value_type> config_ {};

public:
  ConfigParser(const pybind11::handle&);

  template<typename T>
  [[nodiscard]] std::optional<T> get(const std::string& key) const
  {
    if(not config_.contains(key)) { return {}; }

    auto convert_to_requested = []<typename U>(const U& value) -> std::optional<T> {
      if constexpr(std::convertible_to<U, T>) {
        return static_cast<T>(value);
      } else if constexpr(std::same_as<U, none_tag_type>) {
        return {};
      } else {
        throw std::runtime_error {
          "Accessed value is not convertible to requested type."
        };
      }
    };

    return std::visit(convert_to_requested, config_.at(key));
  }

  /// Like get, but throws better exception when optional from get is empty.
  template<typename T>
  [[nodiscard]] T get_or_throw(const std::string& key) const
  {
    if(const auto opt = get<T>(key)) {
      return opt.value();
    } else {
      throw std::runtime_error {
        std::format("Required configuration value missing: {}", key)
      };
    }
  }

  [[nodiscard]] bool contains(const std::string& key) const
  { return this->config_.contains(key); }
};

template<std::size_t halo_size>
struct with_halo_type {};

template<std::size_t halo_size>
static constexpr auto with_halo = with_halo_type<halo_size> {};

/// Gets given parameter from the config parses.
///
/// Throws if any of the numbers are non-positive
/// or if the list is not of expected_length.
template<std::size_t halo_size = 0uz>
std::vector<std::ptrdiff_t>
  get_extent_list(
    const toolbox::ConfigParser& p,
    const std::string& name,
    const std::size_t expected_length,
    with_halo_type<halo_size> = {})
{
  auto x = p.get<std::vector<std::ptrdiff_t>>(name);
  if(not x) {
    throw std::runtime_error { std::format("Config does not contain: {}", name) };
  }

  auto v = std::move(x).value();
  if(v.size() != expected_length) {
    throw std::runtime_error { std::format(
      "Config parameter {} is list of length {} which is not the expected length {}.",
      name,
      v.size(),
      expected_length) };
  }

  for(const auto val: v) {
    if(val <= 0) {
      throw std::runtime_error { std::format(
        "{} is expected to only contain positive integers ({} found)",
        name,
        val) };
    }
  }

  if constexpr(halo_size != 0uz) {
    for(auto& val: v) { val += static_cast<std::ptrdiff_t>(2 * halo_size); }
  }

  return v;
}

}  // namespace toolbox
