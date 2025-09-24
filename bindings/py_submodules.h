#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void tools::bind_tools(py::module& m);
void emf::bind_emf(py::module& m);
void pic::bind_pic(py::module& m);
