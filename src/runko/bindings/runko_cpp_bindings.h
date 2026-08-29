// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tools{ void bind_tools(py::module& m); }
namespace emf  { void bind_emf(  py::module& m); }
namespace pic  { void bind_pic(  py::module& m); }
namespace qed  { void bind_qed(  py::module& m); }


#include <pybind11/numpy.h>

#include <array>
#include <vector>

/// a std::vector as a 1D numpy array (copied)
template<typename T>
auto
  to_ndarray(const std::vector<T>& vec)
{
  const auto shape = std::array { vec.size() };
  auto mda         = py::array_t<T, py::array::c_style>(shape);
  auto mut         = mda.template mutable_unchecked<1>();
  for(auto i = 0uz; i < vec.size(); ++i) { mut(i) = vec[i]; }
  return mda;
}
