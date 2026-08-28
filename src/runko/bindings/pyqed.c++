// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
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

namespace {
/// the physics functions, vectorized over numpy arrays, for tests and analysis
void
  bind_physics(py::module& m)
{
  using qed::Vec3;
  m.attr("ALPHA")    = qed::ALPHA;
  m.attr("BQ_G")     = qed::BQ_G;
  m.attr("SYN_GMIN") = qed::syn::GMIN;
  m.attr("BW_XMIN")  = qed::bw::XMIN;
  m.def("b_from_ellB", py::vectorize([](const double ell_B, const double H) { return qed::b_from_ellB(ell_B, H); }));
  m.def("sigma_kn", py::vectorize([](const double x) { return qed::comp::sigma_kn(x); }));
  m.def(
    "rate_compton",
    py::vectorize([](double kx, double ky, double kz, double ux, double uy, double uz) {
      return qed::comp::rate(Vec3(kx, ky, kz), Vec3(ux, uy, uz));
    }));
  m.def(
    "rate_ann",
    py::vectorize([](double ax, double ay, double az, double bx, double by, double bz) {
      return qed::ann::rate(Vec3(ax, ay, az), Vec3(bx, by, bz));
    }));
  m.def(
    "rate_gg",
    py::vectorize([](double ax, double ay, double az, double bx, double by, double bz) {
      return qed::gg::rate(Vec3(ax, ay, az), Vec3(bx, by, bz));
    }));
  m.def(
    "rate_ssa",
    py::vectorize([](double kx, double ky, double kz, double ux, double uy, double uz, double b, double x_min) {
      return qed::ssa::rate(Vec3(kx, ky, kz), Vec3(ux, uy, uz), b, x_min);
    }));
  m.def("ssa_bound", py::vectorize([](double b, double x_min, double xlo, double glo, double ghi) {
          return qed::ssa::bound(b, x_min, xlo, glo, ghi);
        }));
  m.def("G_syn", py::vectorize([](double ce, double cg) { return qed::ssa::G_syn(ce, cg); }));
  m.def("T_syn", py::vectorize([](double chi) { return qed::syn::T(chi); }));
  m.def("g_syn", py::vectorize([](double chi) { return qed::syn::g(chi); }));
  m.def("P_syn", py::vectorize([](double chi, double r) { return qed::syn::P(chi, r); }));
  m.def("r_syn", py::vectorize([](double z, double chi) { return qed::syn::r_of(z, chi); }));
  m.def("rmean_syn", py::vectorize([](double chi, double rmin) { return qed::syn::rmean(chi, rmin); }));
  m.def("xmean_syn", py::vectorize([](double chi, double gam, bool kinetic, double x_min) {
          return qed::syn::xmean(chi, gam, kinetic, x_min);
        }));
  m.def("rate_syn", py::vectorize([](double ux, double uy, double uz, double b, double x_min, bool kinetic) {
          return qed::syn::rate(Vec3(ux, uy, uz), b, x_min, kinetic);
        }));
  m.def("T_bw", py::vectorize([](double chi) { return qed::bw::T(chi); }));
  m.def("P_bw", py::vectorize([](double chi, double r) { return qed::bw::P(chi, r); }));
  m.def("r_bw", py::vectorize([](double z, double chi) { return qed::bw::r_of(z, chi); }));
  m.def("rate_bw", py::vectorize([](double kx, double ky, double kz, double b) {
          return qed::bw::rate(Vec3(kx, ky, kz), b);
        }));
  m.def("facc_syn", py::vectorize([](double gam, double xm, double acc, double kappa) {
          return qed::facc_syn(gam, xm, acc, kappa);
        }));
  m.def("gamma_cooled", py::vectorize([](double gam, double lam, double xm, double b, double N1) {
          return qed::gamma_cooled(gam, lam, xm, b, N1);
        }));
}
}  // namespace

void
  bind_qed(py::module& m_sub)
{
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");
  py::module m_ph = m_sub.def_submodule("physics", "the QED process functions, vectorized");
  bind_physics(m_ph);

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
    .def("tau_T", &qed::Tile<3>::tau_T)
    .def("interact", &qed::Tile<3>::interact)
    .def("counts", &qed::Tile<3>::counts)
    .def("moved", &qed::Tile<3>::moved)
    .def("report", &qed::Tile<3>::report)
    .def(
      "pair_sum",
      [](qed::Tile<3>& t, const std::string& channel, const bool weight) {
        return t.pair_sum(channel, weight);
      },
      py::arg("channel"),
      py::arg("weight") = false)
    .def(
      "reset_momenta",
      [](qed::Tile<3>& t, const std::string& s, const std::vector<std::array<double, 3>>& u) {
        t.reset_momenta(species_of(s), u);
      });
}

}  // namespace qed
