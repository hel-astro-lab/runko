// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/tile.h"

#include "runko/tools/random.h"
#include "thrust/count.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform_scan.h"
#include "tyvi/mdgrid.h"

#include <format>
#include <stdexcept>

namespace qed {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_grid_idx,
  const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf::Tile<D>(tile_grid_idx, conf),
  pic::Tile<D>(tile_grid_idx, conf),
  sampler_ { Box {}, Params {}, {} }
{
  // ---- the roles of the species, from their mass and charge: a massless species holds the
  // photons, a species of unit mass the leptons (electrons by a negative charge, positrons by a
  // positive one), and heavier species (protons, ions) take no part in the reactions
  auto found = std::array { false, false, false };
  for(auto i = 0uz; i < this->number_of_species(); ++i) {
    const auto a = this->particles(i).args();
    if(a.mass == 0.0) {
      if(found[0]) { throw std::runtime_error { "qed::Tile: more than one massless species." }; }
      sp_[0]   = i;
      found[0] = true;
    } else if(a.mass == 1.0 and a.charge != 0.0) {
      const auto r = a.charge < 0.0 ? 1uz : 2uz;
      if(found[r]) { throw std::runtime_error { "qed::Tile: more than one lepton species of a charge sign." }; }
      sp_[r]   = i;
      found[r] = true;
    } else if(a.mass > 1.0) {
      continue;  // proton or ion: not a QED species
    } else {
      throw std::runtime_error {
        std::format("qed::Tile: species {} with mass {} and charge {} is neither a photon (m = 0), a lepton (m = 1, q != 0) nor an ion (m > 1).", i, a.mass, a.charge)
      };
    }
  }
  if(not(found[0] and found[1] and found[2])) {
    throw std::runtime_error {
      "qed::Tile needs a massless species (photons) and lepton species (m = 1) of both charge signs."
    };
  }
  pops_.reserve(3);
  pops_.emplace_back(true, 0.0);
  pops_.emplace_back(false, this->particles(sp_[1]).args().mass);
  pops_.emplace_back(false, this->particles(sp_[2]).args().mass);

  // ---- the box: tile edge in Compton wavelengths, step in light crossings
  const auto cells = conf.get_or_throw<std::vector<std::ptrdiff_t>>("n_cells_per_tile");
  if(cells[0] != cells[1] or cells[1] != cells[2]) {
    throw std::runtime_error { "qed::Tile needs a cubic tile (n_cells_per_tile equal)." };
  }
  const auto cfl = conf.get_or_throw<double>("cfl");
  const auto box = Box { .H_lamC = conf.get_or_throw<double>("qed_H_lamC"),
                         .dt     = cfl / static_cast<double>(cells[0]),
                         .w_lep  = this->particles(sp_[1]).weight(),
                         .ell_B  = conf.get<double>("qed_ell_B").value_or(0.0),
                         .x_min  = conf.get<double>("qed_x_min").value_or(0.0) };
  if(this->particles(sp_[2]).weight() != box.w_lep) {
    throw std::runtime_error { "qed::Tile: electron and positron weights must agree." };
  }

  // ---- the sampler
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
  par.seed       = static_cast<std::uint64_t>(conf.get<std::ptrdiff_t>("qed_seed").value_or(1));
  if(par.nbin < 1 or par.ncell_edge < 1 or par.nsort < 1 or par.smin < 1) {
    throw std::runtime_error { "qed::Tile: nbin, n_reaction_cells, nsort and smin must be >= 1." };
  }
  auto channels = std::vector<std::string> {};
  if(const auto cs = conf.get<std::vector<std::string>>("qed_channels")) {
    channels = cs.value();
  } else if(const auto c1 = conf.get<std::string>("qed_channels")) {
    // a comma-separated list
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
  sampler_ = Sampler { box, par, channels };
}

template<std::size_t D>
TileView
  Tile<D>::view()
{
  auto t = TileView { .con      = { &this->particles(sp_[0]), &this->particles(sp_[1]), &this->particles(sp_[2]) },
                      .pop      = { &pops_[0], &pops_[1], &pops_[2] },
                      .mins     = layout_mins(),
                      .maxs     = layout_maxs(),
                      .id_base  = {},
                      .ids_used = {} };
  return t;
}

template<std::size_t D>
std::array<typename Tile<D>::value_type, 3>
  Tile<D>::layout_mins() const
{
  return { static_cast<value_type>(this->mins[0]),
           static_cast<value_type>(this->mins[1]),
           static_cast<value_type>(this->mins[2]) };
}

template<std::size_t D>
std::array<typename Tile<D>::value_type, 3>
  Tile<D>::layout_maxs() const
{
  return { static_cast<value_type>(this->maxs[0]),
           static_cast<value_type>(this->maxs[1]),
           static_cast<value_type>(this->maxs[2]) };
}

template<std::size_t D>
void
  Tile<D>::sort_particles()
{
  auto t = view();
  sampler_.layout(t);
}

template<std::size_t D>
void
  Tile<D>::interact()
{
  auto t = view();
  // ids for the products of the step: a block per species, returned unused ones are lost
  constexpr auto block = 1uz << 20uz;
  for(auto s = 0uz; s < 3uz; ++s) { t.id_base[s] = this->reserve_ids(sp_[s], block); }
  sampler_.step(t);
  for(auto s = 0uz; s < 3uz; ++s) {
    if(t.ids_used[s] > block) {
      throw std::runtime_error { "qed::Tile::interact: more than 2^20 products of one species in a step." };
    }
  }
}

template<std::size_t D>
bool
  Tile<D>::laid() const
{
  for(const auto& p: pops_) {
    if(not p.laid()) { return false; }
  }
  return true;
}

template<std::size_t D>
void
  Tile<D>::local_communication_postlude(const int mode)
{
  pic::Tile<D>::local_communication_postlude(mode);
  if(static_cast<runko::comm_mode>(mode) == runko::comm_mode::pic_particle) {
    for(auto& p: pops_) { p.invalidate(); }
  }
}

template<std::size_t D>
void
  Tile<D>::set_particles(
    const Species s,
    const std::vector<std::array<double, 3>>& u,
    const std::vector<double>& w)
{
  if(u.size() != w.size()) {
    throw std::runtime_error { "qed::Tile::set_particles: u and w must have the same length." };
  }
  auto& con      = this->particles(species(s));
  con            = pic::ParticleContainer { con.args() };  // empty, same species
  auto rng       = runko::Rng { runko::rng_key(sampler_.par.seed, sampler_.step(), 7, 0), 0 };
  auto states    = std::vector<runko::ParticleState<double>>(u.size());
  const auto lo  = this->mins;
  const auto hi  = this->maxs;
  for(auto i = 0uz; i < u.size(); ++i) {
    states[i] = runko::ParticleState<double> {
      .pos = { lo[0] + (hi[0] - lo[0]) * rng.uniform(),
               lo[1] + (hi[1] - lo[1]) * rng.uniform(),
               lo[2] + (hi[2] - lo[2]) * rng.uniform() },
      .vel = u[i],
      .id  = 0,
      .wgt = w[i]
    };
  }
  this->inject(species(s), std::move(states));
  pop(s).invalidate();
}

template<std::size_t D>
void
  Tile<D>::reset_momenta(const Species s, const std::vector<std::array<double, 3>>& u)
{
  auto& con    = this->particles(species(s));
  const auto n = con.size();
  auto h       = thrust::host_vector<float>(3 * u.size());
  for(auto i = 0uz; i < u.size(); ++i) {
    for(auto q = 0uz; q < 3uz; ++q) { h[3 * i + q] = static_cast<float>(u[i][q]); }
  }
  const auto d       = thrust::device_vector<float>(h);
  const auto ids_mds = con.ids_mds();
  const auto vel_mds = con.vel_mds();
  const auto w       = tyvi::mdgrid_work {};
  // rank of every slot among the alive ones, in slot order
  auto rank     = thrust::device_vector<runko::index_t>(n);
  const auto i0 = thrust::counting_iterator<runko::index_t>(0);
  thrust::transform_exclusive_scan(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    rank.begin(),
    [=](const runko::index_t i) {
      return ids_mds[i][] != runko::dead_prtc_id ? runko::index_t { 1 } : runko::index_t { 0 };
    },
    runko::index_t { 0 },
    thrust::plus<runko::index_t> {});
  const auto nalive = static_cast<std::size_t>(thrust::count_if(
    w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const runko::index_t i) {
      return ids_mds[i][] != runko::dead_prtc_id;
    }));
  if(nalive != u.size()) {
    throw std::runtime_error { "reset_momenta: one momentum per live LP is needed" };
  }
  const auto rk = thrust::raw_pointer_cast(rank.data());
  const auto dp = thrust::raw_pointer_cast(d.data());
  thrust::for_each(w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const runko::index_t i) {
    if(ids_mds[i][] == runko::dead_prtc_id) { return; }
    for(auto q = 0uz; q < 3uz; ++q) { vel_mds[i][q] = dp[3 * rk[i] + q]; }
  });
  w.wait();
  pop(s).touch_all(con);
}

template<std::size_t D>
double
  Tile<D>::pair_sum(const std::string& channel, const bool weight)
{
  auto t = view();
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(channel == channel_name(static_cast<Channel>(i))) {
      return sampler_.pair_sum(static_cast<Channel>(i), t, weight);
    }
  }
  throw std::runtime_error { "pair_sum: unknown channel " + channel };
}

template<std::size_t D>
void
  Tile<D>::inject_pairs_plaw(const double ell, const double p, const double gmin, const double gmax)
{
  auto t = view();
  for(auto s = 0uz; s < 3uz; ++s) { t.id_base[s] = this->reserve_ids(sp_[s], 1uz << 20uz); }
  sampler_.inject_pairs_plaw(t, ell, p, gmin, gmax);
}

template<std::size_t D>
void
  Tile<D>::inject_photons(const double ell, const double theta, const double w, const bool planck)
{
  auto t = view();
  for(auto s = 0uz; s < 3uz; ++s) { t.id_base[s] = this->reserve_ids(sp_[s], 1uz << 20uz); }
  sampler_.inject_photons(t, ell, theta, w, planck);
}

template<std::size_t D>
void
  Tile<D>::heat_diffusive(const double ell)
{
  auto t = view();
  sampler_.heat_diffusive(t, ell);
}

template<std::size_t D>
void
  Tile<D>::escape_photons()
{
  auto t = view();
  sampler_.escape_photons(t);
}

template<std::size_t D>
double
  Tile<D>::tau_ssa(const double x)
{
  auto t = view();
  return sampler_.tau_ssa(t, x);
}

template<std::size_t D>
std::vector<std::int32_t>
  Tile<D>::get_slot_keys(const Species s)
{
  return pop(s).get_slot_keys(this->particles(species(s)));
}

template<std::size_t D>
Pop::Tables
  Tile<D>::get_tables(const Species s)
{
  return pop(s).get_tables();
}

template<std::size_t D>
Pop::Totals
  Tile<D>::totals(const Species s) const
{
  return Pop::totals(this->particles(species(s)));
}

template<std::size_t D>
double
  Tile<D>::tau_T() const
{
  return totals(Species::em).W + totals(Species::ep).W;
}

template class Tile<3>;

}  // namespace qed
