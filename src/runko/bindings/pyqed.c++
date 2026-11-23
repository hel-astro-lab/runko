// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/actions/qed.h"
#include "runko/bindings/runko_cpp_bindings.h"
#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"

#include <string>
#include <vector>

namespace {

/// the log edges of the escaping spectrum and the linear edges of its direction cosine
std::vector<double>
  esc_edges()
{
  auto e = std::vector<double>(qed::Sampler::NESC + 1);
  for(auto i = 0uz; i < e.size(); ++i) {
    e[i] = qed::Sampler::ESC_LO * std::pow(
                                    qed::Sampler::ESC_HI / qed::Sampler::ESC_LO,
                                    static_cast<double>(i) / qed::Sampler::NESC);
  }
  return e;
}
std::vector<double>
  mu_edges()
{
  auto e = std::vector<double>(qed::Sampler::NMU + 1);
  for(auto i = 0uz; i < e.size(); ++i) {
    e[i] = -1.0 + 2.0 * static_cast<double>(i) / qed::Sampler::NMU;
  }
  return e;
}

/// the species role of its python name
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
  // ---- constants and thresholds
  m.attr("ALPHA")    = qed::ALPHA;
  m.attr("BQ_G")     = qed::BQ_G;
  m.attr("SYN_GMIN") = qed::syn::GMIN;
  m.attr("BW_XMIN")  = qed::bw::XMIN;
  // ---- two-body: cross sections and rates of (incident, partner) momenta
  m.def("b_from_ellB", py::vectorize([](const double ell_B, const double H) {
          return qed::b_from_ellB(ell_B, H);
        }));
  m.def("sigma_kn", py::vectorize([](const double x) {
          return qed::comp::sigma_kn(x);
        }));
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
    py::vectorize([](
                    double kx,
                    double ky,
                    double kz,
                    double ux,
                    double uy,
                    double uz,
                    double b,
                    double x_min) {
      return qed::ssa::rate(Vec3(kx, ky, kz), Vec3(ux, uy, uz), b, x_min);
    }));
  m.def(
    "ssa_bound",
    py::vectorize([](double b, double x_min, double xlo, double glo, double ghi) {
      return qed::ssa::bound(b, x_min, xlo, glo, ghi);
    }));
  // ---- synchrotron: spectral function, rate, cumulative and its inverse, moments
  m.def("G_syn", py::vectorize([](double ce, double cg) {
          return qed::ssa::G_syn(ce, cg);
        }));
  m.def("T_syn", py::vectorize([](double chi) { return qed::syn::T(chi); }));
  m.def("g_syn", py::vectorize([](double chi) { return qed::syn::g(chi); }));
  m.def("P_syn", py::vectorize([](double chi, double r) {
          return qed::syn::P(chi, r);
        }));
  m.def("r_syn", py::vectorize([](double z, double chi) {
          return qed::syn::r_of(z, chi);
        }));
  m.def("rmean_syn", py::vectorize([](double chi, double rmin) {
          return qed::syn::rmean(chi, rmin);
        }));
  m.def(
    "xmean_syn",
    py::vectorize([](double chi, double gam, bool kinetic, double x_min) {
      return qed::syn::xmean(chi, gam, kinetic, x_min);
    }));
  m.def(
    "rate_syn",
    py::vectorize(
      [](double ux, double uy, double uz, double b, double x_min, bool kinetic) {
        return qed::syn::rate(Vec3(ux, uy, uz), b, x_min, kinetic);
      }));
  // ---- Breit-Wheeler: rate, cumulative and its inverse
  m.def("T_bw", py::vectorize([](double chi) { return qed::bw::T(chi); }));
  m.def("P_bw", py::vectorize([](double chi, double r) { return qed::bw::P(chi, r); }));
  m.def("r_bw", py::vectorize([](double z, double chi) {
          return qed::bw::r_of(z, chi);
        }));
  m.def("rate_bw", py::vectorize([](double kx, double ky, double kz, double b) {
          return qed::bw::rate(Vec3(kx, ky, kz), b);
        }));
  // ---- the one-body bookkeeping of the sampler: accumulation factor and cooling law
  m.def("facc_syn", py::vectorize([](double gam, double xm, double acc, double kappa) {
          return qed::facc_syn(gam, xm, acc, kappa);
        }));
  m.def(
    "gamma_cooled",
    py::vectorize([](double gam, double lam, double xm, double b, double N1) {
      return qed::gamma_cooled(gam, lam, xm, b, N1);
    }));
}
}  // namespace

void
  bind_qed(py::module& m_sub)
{
  // runko.qed.threeD.Tile is the python factory of a single QED tile
  // (runko/__init__.py)
  m_sub.def_submodule("threeD", "3D specializations");
  py::module m_ph =
    m_sub.def_submodule("physics", "the QED process functions, vectorized");
  bind_physics(m_ph);
  m_sub.attr("LAMC_CM") = qed::LAMC_CM;
  m_sub.attr("ESC_EDGES") =
    to_ndarray(esc_edges());  // log bins of the escaping spectrum
  m_sub.attr("MU_EDGES") =
    to_ndarray(mu_edges());  // linear bins of its direction cosine

  // ---- one tile of a simulation context, by its id: the per-tile methods of
  // runko.ProxyTile
  using sim_t        = runko::simulation_context;
  using id_t         = qed::tile_id;
  const auto sampler = [](sim_t& sim, const id_t id) -> qed::Sampler& {
    return qed::state(sim, id).sampler;
  };

  // layout: its state, the species roles and the box clocks
  m_sub.def("laid", &qed::laid);
  m_sub.def("layout", &qed::layout);
  m_sub.def("step", [=](sim_t& sim, const id_t id) { return sampler(sim, id).step(); });
  m_sub.def("species", [](sim_t& sim, const id_t id, const std::string& s) {
    return qed::state(sim, id).sp[static_cast<std::size_t>(species_of(s))];
  });
  m_sub.def("box", [=](sim_t& sim, const id_t id) {
    const auto& b = sampler(sim, id).box;
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
  });
  // population IO (host round trips): replace, inspect the layout tables and keys
  m_sub.def(
    "set_particles",
    [](
      sim_t& sim,
      const id_t id,
      const std::string& s,
      const std::vector<std::array<double, 3>>& u,
      const std::vector<double>& w) {
      qed::set_particles(sim, id, species_of(s), u, w);
    });
  m_sub.def("get_slot_keys", [](sim_t& sim, const id_t id, const std::string& s) {
    return to_ndarray(qed::get_slot_keys(sim, id, species_of(s)));
  });
  m_sub.def("get_tables", [](sim_t& sim, const id_t id, const std::string& s) {
    const auto tb =
      qed::state(sim, id).pops[static_cast<std::size_t>(species_of(s))].get_tables();
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
  });
  // diagnostics: totals, pair depth, the sampler step and its counters
  m_sub.def("totals", [](sim_t& sim, const id_t id, const std::string& s) {
    const auto tt = qed::totals(sim, id, species_of(s));
    py::dict d;
    d["N"] = tt.N;
    d["W"] = tt.W;
    d["E"] = tt.E;
    d["P"] = tt.P;
    return d;
  });
  m_sub.def("tau_T", &qed::tau_T);
  m_sub.def("interact", &qed::interact);
  m_sub.def("counts", [=](sim_t& sim, const id_t id) {
    return sampler(sim, id).counts();
  });
  m_sub.def("moved", [=](sim_t& sim, const id_t id) {
    return sampler(sim, id).moved();
  });
  m_sub.def("report", [=](sim_t& sim, const id_t id) {
    return sampler(sim, id).report();
  });
  m_sub.def(
    "pair_sum",
    &qed::pair_sum,
    py::arg("sim"),
    py::arg("id"),
    py::arg("channel"),
    py::arg("weight") = false);
  // the open box: sources, the sink, its budgets and the escaping spectrum
  m_sub.def(
    "inject_pairs_plaw",
    &qed::inject_pairs_plaw,
    py::arg("sim"),
    py::arg("id"),
    py::arg("ell"),
    py::arg("p"),
    py::arg("gmin"),
    py::arg("gmax"));
  m_sub.def(
    "inject_photons",
    &qed::inject_photons,
    py::arg("sim"),
    py::arg("id"),
    py::arg("ell"),
    py::arg("theta"),
    py::arg("w"),
    py::arg("planck") = true);
  m_sub.def(
    "heat_diffusive",
    &qed::heat_diffusive,
    py::arg("sim"),
    py::arg("id"),
    py::arg("ell"));
  m_sub.def(
    "escape_photons",
    &qed::escape_photons,
    py::arg("sim"),
    py::arg("id"),
    py::arg("tau_T")  = -1.0,
    py::arg("dt_esc") = -1.0);
  m_sub.def("escaped_mu_spectrum", [=](sim_t& sim, const id_t id) {
    // (log x bins, mu bins) as a (NESC, NMU) array of x ell_x sums
    auto a = to_ndarray(sampler(sim, id).escaped_mu_spectrum());
    return a.reshape(
      { static_cast<py::ssize_t>(qed::Sampler::NESC),
        static_cast<py::ssize_t>(qed::Sampler::NMU) });
  });
  m_sub.def("reset_budget", [=](sim_t& sim, const id_t id) {
    sampler(sim, id).reset_budget();
  });
  m_sub.def("tau_ssa", &qed::tau_ssa, py::arg("sim"), py::arg("id"), py::arg("x"));
  m_sub.def("budget", [=](sim_t& sim, const id_t id) {
    const auto& b = sampler(sim, id).budget;
    py::dict d;
    d["inj_E"]  = b.inj_E;
    d["inj_W"]  = b.inj_W;
    d["heat_E"] = b.heat_E;
    d["esc_E"]  = b.esc_E;
    d["esc_W"]  = b.esc_W;
    return d;
  });
  m_sub.def("escaped_spectrum", [=](sim_t& sim, const id_t id) {
    // (bin edges, x ell_x per bin) -- the sums are in energy units, the caller divides
    // by steps
    return py::make_tuple(
      to_ndarray(esc_edges()),
      to_ndarray(sampler(sim, id).escaped_spectrum()));
  });
  // in-place momentum reset of the live LPs (thermostats of the tests)
  m_sub.def(
    "reset_momenta",
    [](
      sim_t& sim,
      const id_t id,
      const std::string& s,
      const std::vector<std::array<double, 3>>& u) {
      qed::reset_momenta(sim, id, species_of(s), u);
    });
}

}  // namespace qed
