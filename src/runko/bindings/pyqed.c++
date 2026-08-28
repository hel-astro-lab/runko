// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/tile.h"
#include "runko/qed/virtual_tile.h"
#include "runko/tools/config_parser.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <memory>
#include <vector>

namespace {
namespace py = pybind11;

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

qed::Species
  species_of(const std::string& name)
{
  if(name == "ph") { return qed::Species::ph; }
  if(name == "e-") { return qed::Species::em; }
  if(name == "e+") { return qed::Species::ep; }
  throw std::runtime_error { "qed species must be one of 'ph', 'e-', 'e+'" };
}
}  // namespace

namespace qed {
namespace py = pybind11;

void
  bind_qed(py::module& m_sub)
{
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  py::class_<
    qed::VirtualTile<3>,
    pic::VirtualTile<3>,
    emf::VirtualTile<3>,
    corgi::Tile<3>,
    std::shared_ptr<qed::VirtualTile<3>>>(m_3d, "VirtualTile")
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return qed::VirtualTile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def_static("canonical_type", []() { return py::type::of<qed::Tile<3>>(); });

  py::class_<
    qed::Tile<3>,
    pic::Tile<3>,
    emf::Tile<3>,
    corgi::Tile<3>,
    std::shared_ptr<qed::Tile<3>>>(m_3d, "Tile")
    .def_static("canonical_type", []() { return py::type::of<qed::Tile<3>>(); })
    .def_static(
      "virtual_tile_specialization",
      []() { return py::type::of<qed::VirtualTile<3>>(); })
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return qed::Tile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def("sort_particles", &qed::Tile<3>::sort_particles)
    .def("laid", &qed::Tile<3>::laid)
    .def("step", &qed::Tile<3>::step)
    .def(
      "species",
      [](const qed::Tile<3>& t, const std::string& s) { return t.species(species_of(s)); })
    .def(
      "box",
      [](const qed::Tile<3>& t) {
        const auto& b = t.box();
        py::dict d;
        d["H_lamC"] = b.H_lamC;
        d["dt"]     = b.dt;
        d["w_lep"]  = b.w_lep;
        d["ell_B"]  = b.ell_B;
        d["x_min"]  = b.x_min;
        d["N1"]     = b.N1();
        d["N2"]     = b.N2();
        d["b_perp"] = b.b_perp();
        return d;
      })
    .def(
      "set_particles",
      [](
        qed::Tile<3>& t,
        const std::string& s,
        const std::vector<std::array<double, 3>>& u,
        const std::vector<double>& w) { t.set_particles(species_of(s), u, w); })
    .def(
      "get_slot_keys",
      [](qed::Tile<3>& t, const std::string& s) {
        return to_ndarray(t.get_slot_keys(species_of(s)));
      })
    .def(
      "get_tables",
      [](qed::Tile<3>& t, const std::string& s) {
        const auto tb = t.get_tables(species_of(s));
        py::dict d;
        d["b"]    = to_ndarray(tb.b);
        d["cap"]  = to_ndarray(tb.cap);
        d["fill"] = to_ndarray(tb.fill);
        d["prod"] = to_ndarray(tb.prod);
        d["wmax"] = to_ndarray(tb.wmax);
        d["emin"] = to_ndarray(tb.emin);
        d["emax"] = to_ndarray(tb.emax);
        d["llo"]  = tb.llo;
        d["dl"]   = tb.dl;
        return d;
      })
    .def(
      "totals",
      [](const qed::Tile<3>& t, const std::string& s) {
        const auto tt = t.totals(species_of(s));
        py::dict d;
        d["N"] = tt.N;
        d["W"] = tt.W;
        d["E"] = tt.E;
        d["P"] = tt.P;
        return d;
      })
    .def("tau_T", &qed::Tile<3>::tau_T);
}

}  // namespace qed
