// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tools{ void bind_tools(py::module& m); }
namespace emf  { void bind_emf(  py::module& m); }
namespace pic  { void bind_pic(  py::module& m); }

