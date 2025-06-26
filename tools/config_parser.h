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
};


}  // namespace toolbox
