// Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/particles_common.h"
#include "runko/pic/reflector_wall.h"
#include "runko/pic/tile.h"
#include "runko/pic/virtual_tile.h"
#include "runko/io/pic_average_kinetic_energy.h"
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

namespace {
namespace py = pybind11;
template<typename T>
auto
  to_ndarray(const std::vector<T>& vec)
{
  const auto grid_shape = std::array { vec.size() };
  auto mda              = py::array_t<T, py::array::c_style>(grid_shape);
  auto mda_mut          = mda.template mutable_unchecked<1>();

  for(const auto i: std::views::iota(0uz, vec.size())) { mda_mut(i) = vec[i]; }

  return std::move(mda);
}
}  // namespace

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
    .def(py::init<PS::vec3, PS::vec3>(), py::arg("pos"), py::arg("vel"))
    .def_readwrite("pos", &PS::pos)
    .def_readwrite("vel", &PS::vel);

  // object for storing a batch of multiple particle's data
  using PSB = pic::ParticleStateBatch;
  py::class_<PSB>(m_3d, "ParticleStateBatch")
    .def(
      py::init<PSB::container_type, PSB::container_type>(),
      py::arg("pos"),
      py::arg("vel"))
    .def_readwrite("pos", &PSB::pos)
    .def_readwrite("vel", &PSB::vel);

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

  m_3d.def("_write_average_kinetic_energy", &pic::write_average_kinetic_energy);
}

}  // namespace pic
