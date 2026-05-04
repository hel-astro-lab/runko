// Copyright 2025 - 2026, Joonas Nättilä and the hel-astro-lab contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tools{ void bind_tools(py::module& m); }
namespace emf  { void bind_emf(  py::module& m); }
namespace pic  { void bind_pic(  py::module& m); }

