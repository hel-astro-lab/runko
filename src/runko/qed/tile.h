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
  Sampler sampler_;

  [[nodiscard]] TileView view();

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
  [[nodiscard]] const Params& params() const { return sampler_.par; }
  [[nodiscard]] const Sampler& sampler() const { return sampler_; }
  [[nodiscard]] std::size_t species(const Species s) const
  {
    return sp_[static_cast<std::size_t>(s)];
  }
  [[nodiscard]] Pop& pop(const Species s) { return pops_[static_cast<std::size_t>(s)]; }
  [[nodiscard]] std::size_t step() const { return sampler_.step(); }

  /// The padded reaction layout of every species: the sort of this tile.
  void sort_particles();

  /// One step of the QED Monte Carlo.
  void interact();

  /// Whether every species holds a valid layout.
  [[nodiscard]] bool laid() const;

  /// Particles arrived from other tiles: the layout is stale.
  void local_communication_postlude(const int) override;

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
  [[nodiscard]] std::vector<std::int32_t> get_slot_keys(Species s);
  [[nodiscard]] Pop::Tables get_tables(Species s);
  [[nodiscard]] Pop::Totals totals(Species s) const;
  [[nodiscard]] double tau_T() const;
  [[nodiscard]] double pair_sum(const std::string& channel, bool weight);
  [[nodiscard]] std::map<std::string, double> counts() const { return sampler_.counts(); }
  [[nodiscard]] std::map<std::string, double> moved() const { return sampler_.moved(); }
  [[nodiscard]] std::string report() const { return sampler_.report(); }

  // ---- the open box
  void inject_pairs_plaw(double ell, double p, double gmin, double gmax);
  void inject_photons(double ell, double theta, double w, bool planck);
  void heat_diffusive(double ell);
  void escape_photons();
  [[nodiscard]] const Sampler::Budget& budget() const { return sampler_.budget; }
  void reset_budget() { sampler_.reset_budget(); }
  [[nodiscard]] std::vector<double> escaped_spectrum() const { return sampler_.escaped_spectrum(); }
  [[nodiscard]] double tau_ssa(double x);
};

}  // namespace qed
