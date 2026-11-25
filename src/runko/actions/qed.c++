// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/qed.h"

#include "runko/actions/emf.h"
#include "runko/actions/pic.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/coords.h"
#include "runko/emf/yee_lattice.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/random.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/logical.h"
#include "thrust/sort.h"
#include "thrust/transform_reduce.h"
#include "thrust/transform_scan.h"
#include "tyvi/mdgrid.h"

#include <format>
#include <stdexcept>

namespace qed {

namespace {

template<typename T>
T&
  component(runko::simulation_context& sim, const tile_id id, const char* what)
{
  if(const auto p = sim.tiles.try_get<T>(id)) { return *p; }
  throw std::runtime_error { std::format("qed: the tile has no {}.", what) };
}

pic::particle_containers&
  containers(runko::simulation_context& sim, const tile_id id)
{ return component<pic::particle_containers>(sim, id, "particle containers"); }

runko::global_coordinates_closure<3>
  coordinates(runko::simulation_context& sim, const tile_id id)
{
  const auto& idx = component<runko::cartesian_index<3>>(sim, id, "cartesian index");
  return runko::global_coordinates(sim, idx.as<double>().data);
}

/// Builds the state of a tile from the config: species roles, box and sampler knobs.
tile_state
  make_state(
    const toolbox::ConfigParser& conf,
    const pic::particle_containers& pc,
    const runko::cartesian_index<3>& idx)
{
  // ---- species roles: massless = photons, unit mass = leptons (the sign of the charge
  // picks e- or e+); heavier species (protons, ions) take no part
  auto sp    = std::array<std::size_t, 3> {};
  auto found = std::array { false, false, false };

  for(const auto& [i, con]: pc) {
    const auto a = con.args();
    if(a.mass == 0.0) {
      if(found[0]) {
        throw std::runtime_error { "qed tile: more than one photon species." };
      }
      sp[0]    = i;
      found[0] = true;
    } else if(a.mass == 1.0 and a.charge != 0.0) {
      const auto r = a.charge < 0.0 ? 1uz : 2uz;
      if(found[r]) {
        throw std::runtime_error {
          "qed tile: more than one lepton species of a charge sign."
        };
      }
      sp[r]    = i;
      found[r] = true;
    } else if(a.mass > 1.0) {
      continue;
    } else {
      throw std::runtime_error { std::format(
        "qed tile: species {} with mass {} and charge {} is neither a photon (m = 0), "
        "a lepton (m = 1, q != 0) nor an ion (m > 1).",
        i,
        a.mass,
        a.charge) };
    }
  }

  if(not(found[0] and found[1] and found[2])) {
    throw std::runtime_error {
      "qed tile needs a massless species (photons) and lepton species (m = 1) of both "
      "charge signs."
    };
  }

  if(pc.at(sp[0]).n_payload() < 1uz) {
    throw std::runtime_error { std::format(
      "qed tile: the photons carry their weights as payload: set n_payload{}: 1 for "
      "the photon species.",
      sp[0]) };
  }
  if(pc.at(sp[1]).n_payload() != 0uz or pc.at(sp[2]).n_payload() != 0uz) {
    throw std::runtime_error {
      "qed tile: the leptons carry the species weight w{i}, not a payload."
    };
  }

  // ---- the box: tile edge in Compton wavelengths, step in light crossings
  const auto cells = conf.get_or_throw<std::vector<std::ptrdiff_t>>("n_cells_per_tile");

  // the weight unit W_0 = H^2/sigma_T assumes V = H^3
  if(cells[0] != cells[1] or cells[1] != cells[2]) {
    throw std::runtime_error {
      "qed tile needs a cubic tile (n_cells_per_tile equal)."
    };
  }

  const auto cfl    = conf.get_or_throw<double>("cfl");  // c dt/dx
  const auto tiles  = conf.get_or_throw<std::vector<std::ptrdiff_t>>("n_tiles");
  const auto H_tile = conf.get_or_throw<double>("qed_H_lamC");

  // the compactness is defined on the box, not the tile
  const auto box = Box { .H_lamC     = H_tile,
                         .H_box_lamC = H_tile * static_cast<double>(tiles[0]),
                         .dt         = cfl / static_cast<double>(cells[0]),
                         .w_lep      = runko::species_weight(conf, sp[1]),
                         .ell_B      = conf.get<double>("qed_ell_B").value_or(0.0),
                         .x_min      = conf.get<double>("qed_x_min").value_or(0.0) };
  if(runko::species_weight(conf, sp[2]) != box.w_lep) {
    throw std::runtime_error { "qed tile: electron and positron weights must agree." };
  }

  // ---- the sampler knobs
  const auto geti = [&](const char* key, const int def) {
    return static_cast<int>(conf.get<std::ptrdiff_t>(key).value_or(def));
  };

  auto par       = Params {};
  par.nbin       = geti("qed_nbin", par.nbin);
  par.ncell_edge = geti("qed_n_reaction_cells", par.ncell_edge);
  par.nsort      = geti("qed_nsort", par.nsort);
  par.smin       = geti("qed_smin", par.smin);
  par.n_cap      = geti("qed_n_cap", par.n_cap);
  par.nsub_max   = geti("qed_nsub_max", par.nsub_max);
  par.slack      = conf.get<double>("qed_slack").value_or(par.slack);
  par.a_res      = conf.get<double>("qed_a_res").value_or(par.a_res);
  par.x_ref      = conf.get<double>("qed_x_ref").value_or(par.x_ref);
  par.acc_syn    = conf.get<double>("qed_acc_syn").value_or(par.acc_syn);
  par.kappa      = conf.get<double>("qed_kappa_acc").value_or(par.kappa);
  par.seed =
    static_cast<std::uint64_t>(conf.get<std::ptrdiff_t>("qed_seed").value_or(1));

  // the streams are keyed by seed, step and lane only: fold the grid index into the
  // seed so that no two tiles share a stream
  for(auto q = 0uz; q < 3uz; ++q) {
    par.seed = runko::mix64(par.seed ^ static_cast<std::uint64_t>(idx[q]));
  }

  if(par.nbin < 1 or par.ncell_edge < 1 or par.nsort < 1 or par.smin < 1) {
    throw std::runtime_error {
      "qed tile: nbin, n_reaction_cells, nsort and smin must be >= 1."
    };
  }

  // channels: a list, a comma-separated string, or absent (all the field allows)
  auto channels = std::vector<std::string> {};
  if(const auto cs = conf.get<std::vector<std::string>>("qed_channels")) {
    channels = cs.value();
  } else if(const auto c1 = conf.get<std::string>("qed_channels")) {
    auto rest = c1.value();
    while(not rest.empty()) {
      const auto pos = rest.find(',');
      auto item      = rest.substr(0, pos);
      rest           = pos == std::string::npos ? "" : rest.substr(pos + 1);
      while(not item.empty() and item.front() == ' ') { item.erase(0, 1); }
      while(not item.empty() and item.back() == ' ') { item.pop_back(); }
      if(not item.empty()) { channels.push_back(item); }
    }
  }

  // one layout per role; only the photon keys carry the sub-range bit
  return tile_state { .sp   = sp,
                      .layouts = { Layout { true, 0.0, runko::species_weight(conf, sp[0]) },
                                Layout { false, pc.at(sp[1]).args().mass, box.w_lep },
                                Layout { false, pc.at(sp[2]).args().mass, box.w_lep } },
                      .sampler = Sampler { box, par, channels },
                      .B0      = conf.get<double>("qed_B0").value_or(0.0) };
}

/// Summarizes the lattice field: median |B| over the cells and mean direction.
void
  summarize_field(tile_state& st, emf::YeeLattice& yee)
{
  using I         = runko::index_t;
  const auto Bmds = std::get<1>(yee.view_EBJ_on_device());
  const auto n1 = static_cast<I>(Bmds.extent(1)), n2 = static_cast<I>(Bmds.extent(2));
  const auto n = static_cast<I>(Bmds.extent(0)) * n1 * n2;

  // B of cell f in flat order
  const auto at = [=](const I f) {
    const auto i = f / (n1 * n2), j = (f / n2) % n1, k = f % n2;
    return Vec3(
      static_cast<double>(Bmds[i, j, k][0]),
      static_cast<double>(Bmds[i, j, k][1]),
      static_cast<double>(Bmds[i, j, k][2]));
  };

  const auto w  = tyvi::mdgrid_work {};
  const auto f0 = thrust::counting_iterator<I>(0);
  st.babs.resize(static_cast<std::size_t>(n));
  const auto babs = thrust::raw_pointer_cast(st.babs.data());
  thrust::for_each(w.on_this(), f0, f0 + n, [=](const I f) {
    babs[f] = static_cast<float>(toolbox::norm(at(f)));
  });

  const auto bsum = thrust::transform_reduce(
    w.on_this(),
    f0,
    f0 + n,
    at,
    Vec3(0.0, 0.0, 0.0),
    [](const Vec3& a, const Vec3& b) { return a + b; });
  w.wait();

  thrust::sort(st.babs.begin(), st.babs.end());
  st.b_tile     = static_cast<double>(st.babs[st.babs.size() / 2]);
  const auto nb = toolbox::norm(bsum);
  st.bhat       = nb > 0.0 ? bsum / nb : Vec3(0.0, 0.0, 1.0);
  st.summarized = true;
}

/// Whether any cell of the lattice carries a magnetic field.
bool
  lattice_has_field(emf::YeeLattice& yee)
{
  using I         = runko::index_t;
  const auto Bmds = std::get<1>(yee.view_EBJ_on_device());
  const auto n1 = static_cast<I>(Bmds.extent(1)), n2 = static_cast<I>(Bmds.extent(2));
  const auto n  = static_cast<I>(Bmds.extent(0)) * n1 * n2;
  const auto f0 = thrust::counting_iterator<I>(0);

  return thrust::any_of(tyvi::mdgrid_work {}.on_this(), f0, f0 + n, [=](const I f) {
    const auto i = f / (n1 * n2), j = (f / n2) % n1, k = f % n2;
    return Bmds[i, j, k][0] != 0.0f or Bmds[i, j, k][1] != 0.0f or
           Bmds[i, j, k][2] != 0.0f;
  });
}

/// The field the channels see: uniform, or the lattice scaled to the critical field,
/// with the tile summary of the last interact() (of the first call if none ran).
FieldView
  field_view(
    tile_state& st,
    emf::YeeLattice& yee,
    const runko::global_coordinates_closure<3>& gc)
{
  using yee_value_type = emf::YeeLattice::value_type;
  const auto origo_pos =
    std::array { static_cast<yee_value_type>(gc.mins()[0]) - emf::halo_size,
                 static_cast<yee_value_type>(gc.mins()[1]) - emf::halo_size,
                 static_cast<yee_value_type>(gc.mins()[2]) - emf::halo_size };
  const auto at = yee.interpolate_EB_linear_1st(origo_pos);
  const auto b0 = st.sampler.box.b_perp();

  if(st.B0 <= 0.0 or b0 <= 0.0) {
    return FieldView { .uniform = true,
                       .b0      = b0,
                       .bscale  = 0.0,
                       .b_tile  = b0,
                       .bhat    = Vec3(0.0, 0.0, 1.0),
                       .at      = at };
  }

  if(not st.summarized) { summarize_field(st, yee); }
  return FieldView { .uniform = false,
                     .b0      = b0,
                     .bscale  = b0 / st.B0,
                     .b_tile  = b0 / st.B0 * st.b_tile,
                     .bhat    = st.bhat,
                     .at      = at };
}

/// The tile as the sampler takes it: field, containers, layouts and extent.
TileView
  view(runko::simulation_context& sim, const tile_id id)
{
  auto& st      = state(sim, id);
  auto& pc      = containers(sim, id);
  const auto gc = coordinates(sim, id);
  using F       = pic::ParticleContainer::value_type;

  return TileView {
    .field    = field_view(st, component<emf::YeeLattice>(sim, id, "yee lattice"), gc),
    .con      = { &pc.at(st.sp[0]), &pc.at(st.sp[1]), &pc.at(st.sp[2]) },
    .layout   = { &st.layouts[0], &st.layouts[1], &st.layouts[2] },
    .mins     = runko::cast_array<F>(gc.mins()),
    .maxs     = runko::cast_array<F>(gc.maxs()),
    .id_base  = {},
    .ids_used = {}
  };
}

/// Calls f(sampler, view) with the products numbered from the next free id of their
/// species, then reserves the ids the sampler reports used.
template<typename F>
void
  with_ids(runko::simulation_context& sim, const tile_id id, F&& f)
{
  auto t    = view(sim, id);
  auto& st  = state(sim, id);
  auto& gen = component<pic::particle_id_generator>(sim, id, "particle id generator");

  for(auto s = 0uz; s < 3uz; ++s) { t.id_base[s] = gen.reserve_ids(st.sp[s], 0); }
  f(st.sampler, t);

  for(auto s = 0uz; s < 3uz; ++s) {
    std::ignore = gen.reserve_ids(st.sp[s], t.ids_used[s]);
  }
}

}  // namespace

void
  ensure_constructed_qed(runko::simulation_context& sim)
{
  emf::ensure_constructed_yee_lattices(sim);
  pic::ensure_constructed_particle_containers(sim);

  // collected first: the states are emplaced outside the view that excludes them
  auto todo = std::vector<tile_id> {};
  for(const auto id: sim.tiles.view<pic::particle_containers, runko::local_tile_tag>(
        entt::exclude<tile_state>)) {
    todo.push_back(id);
  }

  for(const auto id: todo) {
    const auto& idx = component<runko::cartesian_index<3>>(sim, id, "cartesian index");
    sim.tiles.emplace<tile_state>(id, make_state(sim.config, containers(sim, id), idx));
  }
}

tile_state&
  state(runko::simulation_context& sim, const tile_id id)
{
  auto& st =
    component<tile_state>(sim, id, "qed state (build it with runko.qed.threeD.Tile)");
  auto& pc = containers(sim, id);

  // the particle exchange appends arrivals behind the last alive slot and so resizes
  // the container; the tables follow it on every access
  for(auto s = 0uz; s < 3uz; ++s) { st.layouts[s].resized(pc.at(st.sp[s])); }
  return st;
}

void
  interact(runko::simulation_context& sim, const tile_id id)
{
  auto& st = state(sim, id);

  // an LP's reaction cell is fixed at the layout and a push between layouts moves it:
  // the step re-widens energy bounds, but nothing repairs a stale cell key
  if(st.sampler.par.ncell_edge > 1 and st.sampler.par.nsort > 1) {
    throw std::runtime_error {
      "qed interact: qed_n_reaction_cells > 1 needs qed_nsort = 1 (a cell key cannot "
      "survive a push)."
    };
  }

  if(st.sampler.box.b_perp() > 0.0) {
    auto& yee = component<emf::YeeLattice>(sim, id, "yee lattice");
    if(st.B0 > 0.0) {
      // the lattice moved since the last lap; only self-absorption reads the summary
      const auto ssa = st.sampler.enabled[static_cast<std::size_t>(Channel::ssa_m)] or
                       st.sampler.enabled[static_cast<std::size_t>(Channel::ssa_p)];
      if(ssa) { summarize_field(st, yee); }
    } else if(lattice_has_field(yee)) {
      throw std::runtime_error {
        "qed interact: the lattice carries a field but qed_B0 is not set; the channels "
        "would see the uniform field of ell_B alone."
      };
    }
  }

  with_ids(sim, id, [](Sampler& s, TileView& t) { s.step(t); });
}

void
  layout(runko::simulation_context& sim, const tile_id id)
{
  auto t = view(sim, id);
  state(sim, id).sampler.layout(t);
}

void
  inject_pairs_plaw(
    runko::simulation_context& sim,
    const tile_id id,
    const double ell,
    const double p,
    const double gmin,
    const double gmax)
{
  with_ids(sim, id, [&](Sampler& s, TileView& t) {
    s.inject_pairs_plaw(t, ell, p, gmin, gmax);
  });
}

void
  inject_photons(
    runko::simulation_context& sim,
    const tile_id id,
    const double ell,
    const double theta,
    const double w,
    const bool planck)
{
  with_ids(sim, id, [&](Sampler& s, TileView& t) {
    s.inject_photons(t, ell, theta, w, planck);
  });
}

void
  heat_diffusive(runko::simulation_context& sim, const tile_id id, const double ell)
{
  auto t = view(sim, id);
  state(sim, id).sampler.heat_diffusive(t, ell);
}

void
  escape_photons(
    runko::simulation_context& sim,
    const tile_id id,
    const double tau_T,
    const double dt_esc)
{
  auto t = view(sim, id);
  state(sim, id).sampler.escape_photons(t, tau_T, dt_esc);
}

void
  set_particles(
    runko::simulation_context& sim,
    const tile_id id,
    const Species s,
    const std::vector<std::array<double, 3>>& u,
    const std::vector<double>& w)
{
  if(u.size() != w.size()) {
    throw std::runtime_error {
      "qed set_particles: u and w must have the same length."
    };
  }

  auto& st      = state(sim, id);
  const auto sp = st.sp[static_cast<std::size_t>(s)];
  auto& con     = containers(sim, id).at(sp);
  auto args     = con.args();
  args.N        = 0;  // args().N is the old size, and fresh slots carry no dead tag
  con           = pic::ParticleContainer { args };

  // the species' own stream: species set in one step must not share positions
  auto rng = runko::Rng { runko::rng_key(
                            st.sampler.par.seed,
                            st.sampler.step(),
                            static_cast<std::uint64_t>(s),
                            stream::positions),
                          0 };

  auto states   = std::vector<runko::ParticleState<double>>(u.size());
  const auto np = con.n_payload();
  auto payload  = std::vector<pic::ParticleContainer::value_type>(u.size() * np, 0.0f);
  const auto gc = coordinates(sim, id);
  const auto lo = gc.mins();
  const auto hi = gc.maxs();
  auto& gen = component<pic::particle_id_generator>(sim, id, "particle id generator");

  for(auto i = 0uz; i < u.size(); ++i) {
    states[i] =
      runko::ParticleState<double> { .pos = { lo[0] + (hi[0] - lo[0]) * rng.uniform(),
                                              lo[1] + (hi[1] - lo[1]) * rng.uniform(),
                                              lo[2] + (hi[2] - lo[2]) * rng.uniform() },
                                     .vel = u[i],
                                     .id  = 0 };
    if(np > 0uz) { payload[i * np] = static_cast<float>(w[i]); }
  }

  // the pic path: ids assigned, then appended
  for(auto& p: states) { p.id = gen.consume_next_id(sp); }
  con.add_particles(states, payload);

  // the container was rebuilt: lay it out again at the next step
  st.layouts[static_cast<std::size_t>(s)].invalidate();
}

void
  reset_momenta(
    runko::simulation_context& sim,
    const tile_id id,
    const Species s,
    const std::vector<std::array<double, 3>>& u)
{
  auto& st     = state(sim, id);
  auto& con    = containers(sim, id).at(st.sp[static_cast<std::size_t>(s)]);
  const auto n = con.size();

  if(con.number_alive() != u.size()) {
    throw std::runtime_error { "reset_momenta: one momentum per live LP is needed" };
  }

  auto h = thrust::host_vector<float>(3 * u.size());
  for(auto i = 0uz; i < u.size(); ++i) {
    for(auto q = 0uz; q < 3uz; ++q) { h[3 * i + q] = static_cast<float>(u[i][q]); }
  }

  const auto d       = thrust::device_vector<float>(h);
  const auto ids_mds = con.ids_mds();
  const auto vel_mds = con.vel_mds();
  const auto w       = tyvi::mdgrid_work {};

  // rank of every slot among the live ones
  auto rank     = thrust::device_vector<runko::index_t>(n);
  const auto i0 = thrust::counting_iterator<runko::index_t>(0);
  thrust::transform_exclusive_scan(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    rank.begin(),
    [=](const runko::index_t i) {
      return ids_mds[i][] != runko::dead_prtc_id ? runko::index_t { 1 }
                                                 : runko::index_t { 0 };
    },
    runko::index_t { 0 },
    thrust::plus<runko::index_t> {});
  const auto rk = thrust::raw_pointer_cast(rank.data());
  const auto dp = thrust::raw_pointer_cast(d.data());

  // a live slot takes the momentum of its rank
  thrust::for_each(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    [=](const runko::index_t i) {
      if(ids_mds[i][] == runko::dead_prtc_id) { return; }
      for(auto q = 0uz; q < 3uz; ++q) { vel_mds[i][q] = dp[3 * rk[i] + q]; }
    });

  // the new energies may lie outside the bounds of their keys
  st.layouts[static_cast<std::size_t>(s)].touch_all(con, w);
  w.wait();
}

bool
  laid(runko::simulation_context& sim, const tile_id id)
{
  auto& st = state(sim, id);
  auto& pc = containers(sim, id);
  for(auto s = 0uz; s < 3uz; ++s) {
    if(not st.layouts[s].laid(pc.at(st.sp[s]))) { return false; }
  }

  return true;
}

std::vector<std::int32_t>
  get_slot_keys(runko::simulation_context& sim, const tile_id id, const Species s)
{
  auto& st = state(sim, id);
  return st.layouts[static_cast<std::size_t>(s)].get_slot_keys(
    containers(sim, id).at(st.sp[static_cast<std::size_t>(s)]));
}

Layout::Totals
  totals(runko::simulation_context& sim, const tile_id id, const Species s)
{
  auto& st = state(sim, id);
  return st.layouts[static_cast<std::size_t>(s)].totals(
    containers(sim, id).at(st.sp[static_cast<std::size_t>(s)]));
}

std::vector<pic::ParticleContainer::value_type>
  get_weights(runko::simulation_context& sim, const tile_id id, const std::size_t ptype)
{
  auto& con = containers(sim, id).at(ptype);
  if(con.n_payload() > 0uz) { return con.get_payload(0); }
  return std::vector(
    con.number_alive(),
    static_cast<pic::ParticleContainer::value_type>(
      runko::species_weight(sim.config, ptype)));
}

double
  tau_T(runko::simulation_context& sim, const tile_id id)
{ return totals(sim, id, Species::em).W + totals(sim, id, Species::ep).W; }

double
  pair_sum(
    runko::simulation_context& sim,
    const tile_id id,
    const std::string& channel,
    const bool weight)
{
  const auto t = view(sim, id);

  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(channel == channel_name(static_cast<Channel>(i))) {
      return state(sim, id).sampler.pair_sum(static_cast<Channel>(i), t, weight);
    }
  }
  throw std::runtime_error { "pair_sum: unknown channel " + channel };
}

double
  tau_ssa(runko::simulation_context& sim, const tile_id id, const double x)
{
  const auto t = view(sim, id);
  return state(sim, id).sampler.tau_ssa(t, x);
}

}  // namespace qed
