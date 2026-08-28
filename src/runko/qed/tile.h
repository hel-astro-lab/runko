// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/pic/tile.h"
#include "runko/qed/box.h"
#include "runko/qed/layout.h"
#include "runko/tools/config_parser.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace qed {

/// The knobs of the sampler.
struct [[nodiscard]] Params {
  int nbin { 16 };          // energy bins per cell and species
  int ncell_edge { 1 };     // reaction sub-cells per tile edge
  int nsort { 10 };         // steps between padded sorts
  int smin { 4 };           // floor on the slack slots of a key
  int n_cap { 32 };         // largest copy number of one reaction
  int nsub_max { 8 };       // sub-steps a channel may take before its strata are let overlap
  double slack { 0.5 };     // slack slots per key, as a fraction of its live LPs
  double a_res { 0.7 };     // exponent of the photon resolution function
  double x_ref { 1e-5 };    // its anchor: one LP per unit weight there
  double acc_syn { 3.0 };   // coefficient of the emission accumulation, 0 disables
  double kappa { 0.1 };     // merged exchange <= kappa (gamma - 1), 0 disables
  std::uint64_t seed { 1 };
};

/// The roles the species play; the container of each is found from its charge and mass.
enum class Species : std::size_t { ph = 0, em = 1, ep = 2 };

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
  Box box_;
  Params par_;
  std::array<std::size_t, 3> sp_ {};  // container index of ph, e-, e+
  std::vector<Pop> pops_;             // ph, e-, e+
  std::size_t istep_ { 0 };

public:
  /// Construct from the config. On top of the pic::Tile keys:
  ///
  /// `qed_H_lamC`: tile edge in Compton wavelengths (the tile must be cubic)
  /// `qed_ell_B`, `qed_x_min`: magnetic compactness and synchrotron floor (0)
  /// `qed_nbin`, `qed_n_reaction_cells`, `qed_nsort`, `qed_slack`, `qed_smin`, `qed_a_res`,
  /// `qed_x_ref`, `qed_n_cap`, `qed_acc_syn`, `qed_kappa_acc`, `qed_nsub_max`, `qed_seed`
  ///
  /// The lepton weight of the box is the electron container's `w{i}`.
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_idx,
    const toolbox::ConfigParser& config);

  ~Tile() = default;

  [[nodiscard]] const Box& box() const { return box_; }
  [[nodiscard]] const Params& params() const { return par_; }
  [[nodiscard]] std::size_t species(const Species s) const
  {
    return sp_[static_cast<std::size_t>(s)];
  }
  [[nodiscard]] Pop& pop(const Species s) { return pops_[static_cast<std::size_t>(s)]; }
  [[nodiscard]] std::size_t step() const { return istep_; }

  /// The padded reaction layout of every species: the sort of this tile.
  void sort_particles();

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
  [[nodiscard]] std::vector<std::int32_t> get_slot_keys(Species s);
  [[nodiscard]] Pop::Tables get_tables(Species s);
  [[nodiscard]] Pop::Totals totals(Species s) const;
  [[nodiscard]] double tau_T() const;
};

}  // namespace qed
