// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

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

/// The QED state of one tile.
///
/// Photons (m = 0) react with electrons (q < 0) and positrons (q > 0) of the tile's
/// particle containers; every tile that takes part carries a tile_state. The reaction
/// layout of a species also sorts its container. The state holds no pointer into the
/// containers (entt relocates components): the sampler's view of the tile is rebuilt
/// on every call.
struct tile_state {
  /// Container index of ph, e-, e+.
  std::array<std::size_t, 3> sp;
  std::array<Layout, 3> layouts;
  Sampler sampler;
  /// Reference field in code units; 0: the uniform field of the box.
  double B0;

  // Tile field of self-absorption, summarized by interact() once per lap: median |B|
  // over the cells (code units) and mean direction.
  double b_tile { 0.0 };
  Vec3 bhat { 0.0, 0.0, 1.0 };
  bool summarized { false };
  /// Scratch of the median.
  thrust::device_vector<float> babs {};
};

/// Construct the yee lattices, particle containers and qed::tile_state of every local
/// tile that lacks them.
///
/// Config keys beyond the pic ones:
///   qed_H_lamC            tile edge in Compton wavelengths (cubic tile)
///   qed_ell_B, qed_x_min  magnetic compactness, synchrotron floor (default 0)
///   qed_channels          channels to run (default: all the field allows)
///   qed_B0                field of ell_B in code units from a PIC driver (0: uniform)
///   qed_seed, qed_nbin, qed_n_reaction_cells, qed_nsort, qed_slack, qed_smin,
///   qed_a_res, qed_x_ref, qed_n_cap, qed_acc_syn, qed_kappa_acc, qed_nsub_max
///
/// The photon species needs `n_payload{i}: 1` (a photon's weight is payload component
/// 0). The leptons carry no payload; both take the electron species' weight `w{i}`.
void ensure_constructed_qed(runko::simulation_context&);

/// Call f(sim, id) on every local tile with a QED state.
///
/// Throws if local tiles exist but none carries one (tiles built as
/// runko.pic.threeD.Tile).
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

// ---- Per tile: the step, the layout, sources and sinks

/// One step of the QED Monte Carlo.
void interact(runko::simulation_context&, tile_id);

/// Lay out every species now; the step otherwise lays out by itself.
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

/// Slab escape at the tile's own depth and step, or, if given, at the box's depth tau_T
/// and step dt_esc (in slab crossings).
void escape_photons(
  runko::simulation_context&,
  tile_id,
  double tau_T  = -1.0,
  double dt_esc = -1.0);

// ---- Per tile: IO and tests

/// Replace the population of a species: momenta u, weights w (ignored for leptons),
/// uniform random positions in the tile.
void set_particles(
  runko::simulation_context&,
  tile_id,
  Species,
  const std::vector<std::array<double, 3>>& u,
  const std::vector<double>& w);

/// Overwrite the momenta of the live LPs of a species in slot order, keeping the
/// layout.
void reset_momenta(
  runko::simulation_context&,
  tile_id,
  Species,
  const std::vector<std::array<double, 3>>& u);

/// Returns the state of a tile, its layouts first resized to the containers.
[[nodiscard]] tile_state& state(runko::simulation_context&, tile_id);

/// Whether every species holds a valid layout.
[[nodiscard]] bool laid(runko::simulation_context&, tile_id);

[[nodiscard]] std::vector<std::int32_t>
  get_slot_keys(runko::simulation_context&, tile_id, Species);

[[nodiscard]] Layout::Totals totals(runko::simulation_context&, tile_id, Species);

/// Weights of the live LPs of container ptype in slot order: payload component 0 if the
/// species has a payload (the photons), else the species weight w{i}.
[[nodiscard]] std::vector<pic::ParticleContainer::value_type>
  get_weights(runko::simulation_context&, tile_id, std::size_t ptype);

/// Pair depth: the summed lepton weight.
[[nodiscard]] double tau_T(runko::simulation_context&, tile_id);

/// Exact expectation of a channel in one step (events, or the weight they move).
[[nodiscard]] double pair_sum(
  runko::simulation_context&,
  tile_id,
  const std::string& channel,
  bool weight);

/// Self-absorption depth of a photon of energy x.
[[nodiscard]] double tau_ssa(runko::simulation_context&, tile_id, double x);

}  // namespace qed
