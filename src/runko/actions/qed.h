// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The QED Monte Carlo on the tiles of a simulation: photons (m = 0) react with
// electrons (q < 0) and positrons (q > 0) of the tile's particle containers. Every tile
// that takes part carries a qed::tile_state; the padded reaction layout of its species
// is also their sort.

#include "runko/comm/cartesian_grid.h"
#include "runko/qed/box.h"
#include "runko/qed/layout.h"
#include "runko/qed/sampler.h"
#include "runko/simulation_context.h"
#include "thrust/device_vector.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace qed {

using tile_id = runko::simulation_context::tile_id_type;

/// The QED state of one tile. It holds no pointer into the containers (entt relocates
/// components): the sampler's view of the tile is rebuilt on every call.
struct tile_state {
  std::array<std::size_t, 3> sp;  // container index of ph, e-, e+
  std::array<Pop, 3> pops;        // their layouts
  Sampler sampler;                // the Monte Carlo state of the tile
  double B0;  // reference field in code units; 0: the uniform field of the box
  // the tile field of self-absorption, summarized once per lap by interact(): the
  // median |B| over the cells in units of B_Q and the mean direction
  double b_tile { 0.0 };
  Vec3 bhat { 0.0, 0.0, 1.0 };
  bool summarized { false };
  thrust::device_vector<float> babs {};  // |B| per cell, scratch of the median
};

/// Construct the yee lattices, the particle containers and the qed::tile_state of every
/// local tile that lacks them. On top of the pic keys the config gives
///
/// `qed_H_lamC`: tile edge in Compton wavelengths (the tile must be cubic)
/// `qed_ell_B`, `qed_x_min`: magnetic compactness and synchrotron floor (0)
/// `qed_nbin`, `qed_n_reaction_cells`, `qed_nsort`, `qed_slack`, `qed_smin`,
/// `qed_a_res`, `qed_x_ref`, `qed_n_cap`, `qed_acc_syn`, `qed_kappa_acc`,
/// `qed_nsub_max`, `qed_seed` `qed_channels`: the channels to run (all that the field
/// allows if absent) `qed_B0`: the field of ell_B in code units, set by a PIC driver
/// (0: the uniform field)
///
/// The lepton weight of the box is the electron container's `w{i}`.
void ensure_constructed_qed(runko::simulation_context&);

/// f(sim, id) on every local tile with a QED state. A QED action on local tiles none of
/// which carries one (tiles built as runko.pic.threeD.Tile) is a driver error, not a
/// no-op.
template<typename F>
void
  for_each_tile(runko::simulation_context& sim, F&& f)
{
  auto any = false;
  for(auto&& [id, _]: sim.view_tiles<tile_state, runko::local_tile_tag>()) {
    f(sim, id);
    any = true;
  }
  if(not any and not sim.tiles.view<runko::local_tile_tag>().empty()) {
    throw std::runtime_error {
      "qed: no local tile carries a QED state (build the tiles with "
      "runko.qed.threeD.Tile)."
    };
  }
}

// ---- one tile: the step, its layout, the sources and sinks of the open box
/// One step of the QED Monte Carlo.
void interact(runko::simulation_context&, tile_id);
/// The padded reaction layout of every species: the sort of a QED tile.
void layout(runko::simulation_context&, tile_id);
void inject_pairs_plaw(
  runko::simulation_context&,
  tile_id,
  double ell,
  double p,
  double gmin,
  double gmax);
void inject_photons(
  runko::simulation_context&,
  tile_id,
  double ell,
  double theta,
  double w,
  bool planck);
void heat_diffusive(runko::simulation_context&, tile_id, double ell);
/// slab escape: at this tile's own depth and step, or at the depth tau_T and the step
/// dt_esc (in crossings of the slab) the driver gives for the whole box
void escape_photons(
  runko::simulation_context&,
  tile_id,
  double tau_T  = -1.0,
  double dt_esc = -1.0);

// ---- one tile: IO and tests
/// Replace the population of a species: momenta u, weights w (ignored for leptons),
/// uniform random positions in the tile.
void set_particles(
  runko::simulation_context&,
  tile_id,
  Species,
  const std::vector<std::array<double, 3>>& u,
  const std::vector<double>& w);
/// Redraw the momenta of the live LPs of a species in place (slot order), keeping the
/// layout.
void reset_momenta(
  runko::simulation_context&,
  tile_id,
  Species,
  const std::vector<std::array<double, 3>>& u);
/// The state of a tile, its arrivals keyed as tail first.
[[nodiscard]] tile_state& state(runko::simulation_context&, tile_id);
/// Whether every species holds a valid layout.
[[nodiscard]] bool laid(runko::simulation_context&, tile_id);
[[nodiscard]] std::vector<std::int32_t>
  get_slot_keys(runko::simulation_context&, tile_id, Species);
[[nodiscard]] Pop::Totals totals(runko::simulation_context&, tile_id, Species);
/// pair depth: the lepton weight
[[nodiscard]] double tau_T(runko::simulation_context&, tile_id);
/// exact expectation of a channel in one step (events, or the weight they move)
[[nodiscard]] double pair_sum(
  runko::simulation_context&,
  tile_id,
  const std::string& channel,
  bool weight);
/// self-absorption depth of a photon of energy x
[[nodiscard]] double tau_ssa(runko::simulation_context&, tile_id, double x);

}  // namespace qed
