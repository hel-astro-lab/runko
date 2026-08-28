// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/tile.h"

#include "runko/tools/random.h"

#include <format>
#include <stdexcept>

namespace qed {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_grid_idx,
  const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf::Tile<D>(tile_grid_idx, conf),
  pic::Tile<D>(tile_grid_idx, conf)
{
  // ---- the roles of the species, from their charge and mass
  auto found = std::array { false, false, false };
  for(auto i = 0uz; i < this->number_of_species(); ++i) {
    const auto a = this->particles(i).args();
    if(a.mass == 0.0) {
      sp_[0]   = i;
      found[0] = true;
    } else if(a.charge < 0.0) {
      sp_[1]   = i;
      found[1] = true;
    } else if(a.charge > 0.0) {
      sp_[2]   = i;
      found[2] = true;
    }
  }
  if(not(found[0] and found[1] and found[2])) {
    throw std::runtime_error {
      "qed::Tile needs a massless species (photons), a negative and a positive lepton species."
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
  box_ = Box { .H_lamC = conf.get_or_throw<double>("qed_H_lamC"),
               .dt     = cfl / static_cast<double>(cells[0]),
               .w_lep  = this->particles(sp_[1]).weight(),
               .ell_B  = conf.get<double>("qed_ell_B").value_or(0.0),
               .x_min  = conf.get<double>("qed_x_min").value_or(0.0) };
  if(this->particles(sp_[2]).weight() != box_.w_lep) {
    throw std::runtime_error { "qed::Tile: electron and positron weights must agree." };
  }

  // ---- the sampler
  const auto geti = [&](const char* key, const int def) {
    return static_cast<int>(conf.get<std::ptrdiff_t>(key).value_or(def));
  };
  par_.nbin       = geti("qed_nbin", par_.nbin);
  par_.ncell_edge = geti("qed_n_reaction_cells", par_.ncell_edge);
  par_.nsort      = geti("qed_nsort", par_.nsort);
  par_.smin       = geti("qed_smin", par_.smin);
  par_.n_cap      = geti("qed_n_cap", par_.n_cap);
  par_.nsub_max   = geti("qed_nsub_max", par_.nsub_max);
  par_.slack      = conf.get<double>("qed_slack").value_or(par_.slack);
  par_.a_res      = conf.get<double>("qed_a_res").value_or(par_.a_res);
  par_.x_ref      = conf.get<double>("qed_x_ref").value_or(par_.x_ref);
  par_.acc_syn    = conf.get<double>("qed_acc_syn").value_or(par_.acc_syn);
  par_.kappa      = conf.get<double>("qed_kappa_acc").value_or(par_.kappa);
  par_.seed       = static_cast<std::uint64_t>(conf.get<std::ptrdiff_t>("qed_seed").value_or(1));
  if(par_.nbin < 1 or par_.ncell_edge < 1 or par_.nsort < 1 or par_.smin < 1) {
    throw std::runtime_error { "qed::Tile: nbin, n_reaction_cells, nsort and smin must be >= 1." };
  }
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
  const auto lp = LayoutParams { .nbin       = par_.nbin,
                                 .ncell_edge = par_.ncell_edge,
                                 .slack      = par_.slack,
                                 .smin       = par_.smin };
  for(auto s = 0uz; s < pops_.size(); ++s) {
    pops_[s].layout(
      this->particles(sp_[s]),
      lp,
      layout_mins(),
      layout_maxs(),
      runko::rng_key(par_.seed, istep_, s, 0));
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
  auto rng       = runko::Rng { runko::rng_key(par_.seed, istep_, 7, 0), 0 };
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
