// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/pic/tile.h"
#include "runko/qed/box.h"
#include "runko/qed/layout.h"
#include "runko/qed/sampler.h"
#include "runko/tools/config_parser.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace qed {

/*! \brief QED tile
 *
 * A pic tile whose species react through the QED Monte Carlo: photons (m = 0) with
 * electrons (q < 0) and positrons (q > 0). The padded reaction layout of every species is its
 * sort, so sort_particles() is the sort the pusher needs and the epoch of the sampler.
 */
template<std::size_t D>
class Tile : virtual public pic::Tile<D> {
public:
  static_assert(D == 3);
  using value_type = pic::ParticleContainer::value_type;

private:
  std::array<std::size_t, 3> sp_ {};  // container index of ph, e-, e+
  std::vector<Pop> pops_;             // ph, e-, e+
  Sampler sampler_;                   // the Monte Carlo state of this tile
  double B0_ { 0.0 };                 // reference field in code units; 0: the uniform field of the box
  // the tile field of self-absorption, summarized once per lap by interact(): the median |B|
  // over the cells in units of B_Q and the mean direction
  double b_tile_ { 0.0 };
  Vec3 bhat_ { 0.0, 0.0, 1.0 };
  bool summarized_ { false };
  thrust::device_vector<float> babs_;   // |B| per cell, scratch of the median

  /// the containers, layouts and extent of this tile as the sampler takes them
  [[nodiscard]] TileView view();
  /// the field the channels see: uniform, or the lattice scaled to the critical field with the
  /// tile summary of the last interact() (of the first call if none ran)
  [[nodiscard]] FieldView field_view();
  /// the median |B| and the mean direction of the lattice field
  void summarize_field();
  /// whether any cell of the lattice carries a magnetic field
  [[nodiscard]] bool lattice_has_field();

public:
  /// Construct from the config. On top of the pic::Tile keys:
  ///
  /// `qed_H_lamC`: tile edge in Compton wavelengths (the tile must be cubic)
  /// `qed_ell_B`, `qed_x_min`: magnetic compactness and synchrotron floor (0)
  /// `qed_nbin`, `qed_n_reaction_cells`, `qed_nsort`, `qed_slack`, `qed_smin`, `qed_a_res`,
  /// `qed_x_ref`, `qed_n_cap`, `qed_acc_syn`, `qed_kappa_acc`, `qed_nsub_max`, `qed_seed`
  /// `qed_channels`: the channels to run (all that the field allows if absent)
  ///
  /// The lepton weight of the box is the electron container's `w{i}`.
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_idx,
    const toolbox::ConfigParser& config);

  ~Tile() = default;

  [[nodiscard]] const Box& box() const { return sampler_.box; }
  [[nodiscard]] std::size_t species(const Species s) const
  {
    return sp_[static_cast<std::size_t>(s)];
  }
  [[nodiscard]] Pop& pop(const Species s) { return pops_[static_cast<std::size_t>(s)]; }
  [[nodiscard]] std::size_t step() const { return sampler_.step(); }

  /// The padded reaction layout of every species: the sort of this tile.
  void sort_particles() override;

  /// One step of the QED Monte Carlo.
  void interact();

  /// Whether every species holds a valid layout.
  [[nodiscard]] bool laid() const;

  /// Particles arrived from other tiles land in the tail of the layout, inert until the next
  /// sort, so the layout stays valid.
  void append_arrivals(
    std::size_t ptype,
    const std::vector<std::span<const runko::ParticleState<value_type>>>& spans,
    std::array<value_type, 3> mins,
    std::array<value_type, 3> maxs) override;

  /// tile extent as the layout needs it
  [[nodiscard]] std::array<value_type, 3> layout_mins() const;
  [[nodiscard]] std::array<value_type, 3> layout_maxs() const;

  // ---- IO and tests
  /// Replace the population of a species: momenta u, weights w (ignored for leptons), uniform
  /// random positions in the tile.
  void set_particles(
    Species s,
    const std::vector<std::array<double, 3>>& u,
    const std::vector<double>& w);
  /// Redraw the momenta of the live LPs of a species in place (slot order), keeping the layout.
  void reset_momenta(Species s, const std::vector<std::array<double, 3>>& u);
  [[nodiscard]] std::vector<std::int32_t> get_slot_keys(Species s);  // key of every live LP
  [[nodiscard]] Pop::Tables get_tables(Species s);                   // the per-key tables
  [[nodiscard]] Pop::Totals totals(Species s) const;                 // N, W, E, P of a species
  [[nodiscard]] double tau_T() const;                                // pair depth: the lepton weight
  /// exact expectation of a channel in one step (events, or the weight they move)
  [[nodiscard]] double pair_sum(const std::string& channel, bool weight);
  [[nodiscard]] std::map<std::string, double> counts() const { return sampler_.counts(); }
  [[nodiscard]] std::map<std::string, double> moved() const { return sampler_.moved(); }
  [[nodiscard]] std::string report() const { return sampler_.report(); }

  // ---- the open box: sources and sinks, each accumulating into the budget
  void inject_pairs_plaw(double ell, double p, double gmin, double gmax);  // pairs, dN/dgamma ~ gamma^-p
  void inject_photons(double ell, double theta, double w, bool planck);    // Wien or Planck photons
  void heat_diffusive(double ell);                                        // momentum diffusion of the leptons
  /// slab escape: at this tile's own depth and step, or at the depth tau_T and the step dt_esc
  /// (in crossings of the slab) the driver gives for the whole box
  void escape_photons(double tau_T = -1.0, double dt_esc = -1.0);
  [[nodiscard]] const Sampler::Budget& budget() const { return sampler_.budget; }
  void reset_budget() { sampler_.reset_budget(); }
  [[nodiscard]] std::vector<double> escaped_spectrum() const { return sampler_.escaped_spectrum(); }
  [[nodiscard]] std::vector<double> escaped_mu_spectrum() const { return sampler_.escaped_mu_spectrum(); }
  [[nodiscard]] double tau_ssa(double x);  // self-absorption depth of a photon of energy x
};

}  // namespace qed
