// Copyright 2025 - 2026, Joonas Nättilä and the hel-astro-lab contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>

#include "core/communication_common.h"
#include "runko_cpp_bindings.h"

//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(runko_cpp_bindings, m_base) {
  m_base.doc() = "Runko Python3 bindings";

  // Not used directly but required to be loaded,
  // in order for python to know about corgi base classes.
  std::ignore = py::module::import("pycorgi");

  /// auxiliary tools
  py::module m_tools = m_base.def_submodule("tools", "auxiliary tools");
  tools::bind_tools(m_tools);

  /// emf
  py::module m_emf = m_base.def_submodule("emf", "emf module");
  emf::bind_emf(m_emf);

  /// pic
  py::module m_pic = m_base.def_submodule("pic", "pic module");
  pic::bind_pic(m_pic);
}

