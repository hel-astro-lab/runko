#pragma once

#include "pybind11/pybind11.h"

#include <concepts>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace toolbox {

class [[nodiscard]] ConfigParser {
  using value_type =
    std::variant<bool, std::string, std::ptrdiff_t, double, std::vector<std::string>>;
  std::unordered_map<std::string, value_type> config_ {};

public:
  ConfigParser(const pybind11::handle&);

  template<typename T>
  [[nodiscard]] std::optional<T> get(const std::string& key) const
  {
    if(not config_.contains(key)) { return {}; }

    auto convert_to_requested = []<typename U>(const U& value) -> T {
      if constexpr(std::convertible_to<U, T>) {
        return static_cast<T>(value);
      } else {
        throw std::runtime_error {
          "Accessed value is not convertible to requested type."
        };
      }
    };

    return std::visit(convert_to_requested, config_.at(key));
  }
};


}  // namespace toolbox
