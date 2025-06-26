#include "tools/config_parser.h"

#include "pybind11/stl.h"

#include <format>
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

      // Assume list is list of strings.
      config_.emplace(
        std::make_pair(std::move(skey), value.cast<std::vector<std::string>>()));
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
