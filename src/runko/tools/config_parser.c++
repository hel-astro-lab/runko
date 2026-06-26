// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/tools/config_parser.h"

#include "pybind11/stl.h"

#include <format>
#include <ranges>
#include <utility>

namespace toolbox {

ConfigParser::ConfigParser(const pybind11::handle& conf_obj)
{
  if(not conf_obj) {
    throw std::runtime_error { "Given configuration object is not valid." };
  }

  if(not pybind11::hasattr(conf_obj, "__dict__")) {
    throw std::runtime_error {
      "Given configuration object does not have __dict__ attribute."
    };
  }

  const pybind11::dict d = conf_obj.attr("__dict__");
  for(const auto& [key, value]: d) {
    auto skey = key.cast<std::string>();
    if(pybind11::isinstance<pybind11::none>(value)) {
      config_.emplace(std::make_pair(std::move(skey), none_tag_type {}));
    } else if(pybind11::isinstance<pybind11::bool_>(value)) {
      config_.emplace(std::make_pair(std::move(skey), value.cast<bool>()));
    } else if(pybind11::isinstance<pybind11::str>(value)) {
      config_.emplace(std::make_pair(std::move(skey), value.cast<std::string>()));
    } else if(pybind11::isinstance<pybind11::int_>(value)) {
      config_.emplace(std::make_pair(std::move(skey), value.cast<std::ptrdiff_t>()));
    } else if(pybind11::isinstance<pybind11::float_>(value)) {
      config_.emplace(std::make_pair(std::move(skey), value.cast<double>()));
    } else if(pybind11::isinstance<pybind11::list>(value)) {

      auto contains_ints    = false;
      auto contains_floats  = false;
      auto contains_strings = false;


      for(const auto& [n, x]:
          std::views::enumerate(pybind11::reinterpret_borrow<pybind11::list>(value))) {
        const auto is_int   = pybind11::isinstance<pybind11::int_>(x);
        const auto is_float = pybind11::isinstance<pybind11::float_>(x);
        const auto is_str   = pybind11::isinstance<pybind11::str>(x);

        if(not(is_int or is_float or is_str)) {
          const auto cls      = x.attr("__class__");
          const auto cls_name = cls.attr("__name__").cast<std::string>();
          const auto msg      = std::format(
            "list element of `{}` at index {} has unsupported type: {}",
            skey,
            n,
            cls_name);
          throw std::runtime_error { msg };
        }

        contains_ints    = contains_ints or is_int;
        contains_floats  = contains_floats or is_float;
        contains_strings = contains_strings or is_str;
      }

      if(contains_ints and contains_floats and contains_strings) {
        throw std::runtime_error { std::format(
          "Can not deduce common type for list that contains ints, floats and strings "
          "(config parameter: {})!",
          skey) };
      } else if(not contains_ints and not contains_floats and contains_strings) {
        config_.emplace(
          std::make_pair(std::move(skey), value.cast<std::vector<std::string>>()));
      } else if(contains_ints and not contains_floats and not contains_strings) {
        config_.emplace(
          std::make_pair(std::move(skey), value.cast<std::vector<std::ptrdiff_t>>()));
      } else {
        config_.emplace(
          std::make_pair(std::move(skey), value.cast<std::vector<double>>()));
      }
    } else {
      const auto cls      = value.attr("__class__");
      const auto cls_name = cls.attr("__name__").cast<std::string>();
      const auto msg      = std::format(
        "configuration variable `{}` has unsupported type: {}",
        skey,
        cls_name);
      throw std::runtime_error { msg };
    }
  }
}
}  // namespace toolbox
