#pragma once

#include "pybind11/pybind11.h"

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
  [[nodiscard]] T get(const std::string& key) const
  {
    try {
      return std::get<T>(config_.at(key));
    } catch(std::out_of_range) {
      throw std::runtime_error { std::string { "Trying to access missing key: " } +
                                 key };
    }
  }
};


}  // namespace toolbox
