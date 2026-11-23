// Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/particles_common.h"
#include "runko/pic/reflector_wall.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/tools/config_parser.h"

#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

//--------------------------------------------------
// experimental PIC module


namespace pic {
namespace py = pybind11;

// python bindings for plasma classes & functions
void
  bind_pic(py::module& m_sub)
{
  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  // object for storing single particle data
  using PS = runko::ParticleState<double>;
  py::class_<PS>(m_3d, "ParticleState")
    .def(
      py::init([](const PS::vec3 pos, const PS::vec3 vel, const double wgt) {
        return PS { .pos = pos, .vel = vel, .id = 0, .wgt = wgt };
      }),
      py::arg("pos"),
      py::arg("vel"),
      py::arg("wgt") = 1.0)
    .def_readwrite("pos", &PS::pos)
    .def_readwrite("vel", &PS::vel)
    .def_readwrite("wgt", &PS::wgt);

  // object for storing a batch of multiple particle's data
  using PSB = runko::ParticleStateBatch;
  py::class_<PSB>(m_3d, "ParticleStateBatch")
    .def(
      py::init([](
                 PSB::container_type pos,
                 PSB::container_type vel,
                 std::optional<py::array_t<double>> wgt,
                 py::object payload) {
        return PSB { .pos     = std::move(pos),
                     .vel     = std::move(vel),
                     .wgt     = std::move(wgt),
                     .payload = std::move(payload) };
      }),
      py::arg("pos"),
      py::arg("vel"),
      py::arg("wgt")     = py::none(),
      py::arg("payload") = py::none())
    .def_readwrite("pos", &PSB::pos)
    .def_readwrite("vel", &PSB::vel)
    .def_readwrite("wgt", &PSB::wgt)
    .def_readwrite("payload", &PSB::payload);

  // reflector wall data structure
  using RW = pic::reflector_wall;
  py::class_<RW>(m_3d, "reflector_wall")
    .def(
      py::init(
        [](RW::value_type walloc, RW::value_type betawall, RW::value_type gammawall) {
          return RW { walloc, betawall, gammawall };
        }),
      py::arg("walloc"),
      py::arg("betawall")  = RW::value_type { 0 },
      py::arg("gammawall") = RW::value_type { 1 })
    .def_readwrite("walloc", &RW::walloc)
    .def_readwrite("betawall", &RW::betawall)
    .def_readwrite("gammawall", &RW::gammawall);
}

}  // namespace pic
