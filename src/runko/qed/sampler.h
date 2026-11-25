// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Cell-envelope Monte Carlo of one tile over the energy-binned padded layout of every
// species. Two-body channels (Compton, pair annihilation, photon-photon pair
// production, synchrotron self-absorption) draw per (cell, bin a, bin b) the envelope
// Lambda = c_l N_a N_b 2 sigmahat_ab W_max/N_2 by stratified pairing (twobody.c++);
// one-body channels (quantum synchrotron emission, magnetic pair production) run one
// lane per slot (onebody.c++). What an event moves and how often it fires is the rule
// of weights.h. Products go into the slack of their key (reactive at once) or the tail
// (inert until the next layout).

#include "runko/pic/particle.h"
#include "runko/qed/box.h"
#include "runko/qed/channels.h"
#include "runko/qed/field.h"
#include "runko/qed/layout.h"
#include "runko/qed/species.h"
#include "runko/qed/weights.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace qed {

/// The knobs of the sampler.
struct [[nodiscard]] Params {
  int nbin { 16 };         // energy bins per cell and species
  int ncell_edge { 1 };    // reaction sub-cells per tile edge
  int nsort { 10 };        // steps between padded sorts
  int smin { 4 };          // floor on the slack slots of a key
  int n_cap { 32 };        // largest copy number of one reaction
  int nsub_max { 8 };      // sub-steps of an additive channel before rows are clamped
  double slack { 0.5 };    // slack slots per key, as a fraction of its live LPs
  double a_res { 0.7 };    // exponent of the photon resolution function
  double x_ref { 1e-5 };   // its anchor: one LP per unit weight there
  double acc_syn { 3.0 };  // coefficient of the emission accumulation, 0 disables
  double kappa { 0.1 };    // merged exchange <= kappa (gamma - 1), 0 disables
  std::uint64_t seed { 1 };
};

/// Per-channel counters of a run.
struct [[nodiscard]] ChannelStats {
  double M { 0 };        // envelope lanes (two-body), slots (one-body)
  double virt { 0 };     // two-body lanes on a dead slot
  double wpass { 0 };    // lanes past the weight stage
  double acc { 0 };      // accepted events
  double wmoved { 0 };   // weight moved, sum f s
  double facc { 0 };     // summed accumulation factors of accepted events
  double wem { 0 };      // emitted weight (one-body emission)
  double clamped { 0 };  // events dropped to keep the strata of a bin distinct
  double shbar_num { 0 }, shbar_den { 0 };  // <sigmahat> of the envelope
  int nsub_max_seen { 1 };                  // largest sub-cycle count of a step
  double nsub_steps { 0 };                  // steps with more than one sub-cycle
};

/// Counters of a run.
struct [[nodiscard]] Stats {
  std::size_t steps { 0 };
  double cell_lanes { 0 };  // (cell, bin, bin) lanes of the envelope kernels
  std::array<ChannelStats, N_CHANNELS> ch {};
  double slackprod { 0 }, tailprod { 0 };   // products placed in slack / in the tail
  double heldprod { 0 };                    // tail ones whose slack held an arrival
  std::array<double, 3> slots {}, live {};  // per species
};

/// A product of a channel pass, staged before placement.
struct [[nodiscard]] Staged {
  float u[3];            // momentum
  float pos[3];          // position (the parent's)
  float w;               // weight of each copy
  std::int32_t n_r;      // copies to place
  std::int32_t cell;     // reaction cell
  std::int32_t species;  // 0 ph, 1 e-, 2 e+
  std::int32_t bit;      // photon sub-range bit
};

/// Appends product records to the staging buffer from inside a kernel.
struct [[nodiscard]] Stager {
  Staged* rec;
  std::int32_t* count;  // records staged so far

  /// Stages n_r copies of weight w of a product of species s at pos, in the reaction
  /// cell `cell` and the photon sub-range `bit`.
  constexpr void operator()(
    const Vec3& u,
    const float* pos,
    const double w,
    const int n_r,
    const int cell,
    const Species s,
    const int bit) const
  {
    const auto k =
      static_cast<std::size_t>(sstd::atomic_fetch_add(count, std::int32_t { 1 }));
    rec[k] = Staged { .u       = { static_cast<float>(u[0]),
                                   static_cast<float>(u[1]),
                                   static_cast<float>(u[2]) },
                      .pos     = { pos[0], pos[1], pos[2] },
                      .w       = static_cast<float>(w),
                      .n_r     = n_r,
                      .cell    = cell,
                      .species = static_cast<std::int32_t>(s),
                      .bit     = bit };
  }
};

/// Sampler state of a tile: scratch buffers, rounding carries and counters kept over
/// steps.
class [[nodiscard]] Sampler {
public:
  Box box;
  Params par;
  std::array<bool, N_CHANNELS> enabled {};  // per Channel
  Stats stats;

private:
  template<typename Ch>
  friend class TwoBody;

  // ---- scratch, kept between passes for its allocations

  /// One side of a two-body channel, as of the start of the channel: per (cell, bin)
  /// the first slot and the fill of its one or two sub-ranges (the second at ncb + cb)
  /// and the bounds over them.
  struct SideTable {
    thrust::device_vector<runko::index_t> base, cnt;
    thrust::device_vector<float> wmax, emin, emax;
  };
  std::array<SideTable, 2> sides_;  // incident, partner

  /// Envelope tables per (cell, bin a, bin b): the envelope, the rate bound and the
  /// weight bound; and the envelope summed over b (row) and over a (col).
  struct Envelope {
    thrust::device_vector<double> Lam, bound, wmax, row, col;
  } env_;

  /// Envelope events of a sub-step per (cell, bin a, bin b): the count M and its
  /// exclusive scan, the first lane; for a consumed partner the strata base along a
  /// and the events per partner bin.
  struct Strata {
    thrust::device_vector<std::int32_t> M, off2, M2;
    thrust::device_vector<runko::index_t> first;
  } strata_;

  /// One envelope event: the slot pair it proposes, the weight its event would move,
  /// its bin pair and how the thinning ended.
  struct Lane {
    runko::index_t i, j;
    double S;
    std::int32_t g;
    std::uint8_t state;
  };
  thrust::device_vector<Lane> lanes_;

  thrust::device_vector<runko::index_t> acc_;  // accepted lanes; fired slots (one-body)
  thrust::device_vector<std::uint8_t> fired_;  // one-body: the slot fired

  /// Accepted events: the probabilities q1, q2 that parent 1, 2 is updated or removed
  /// and the pair count q3, their rounded decisions n1, n2, n3, the cell, the exchanges
  /// f the partner takes as one and its momentum kick (3 per event).
  struct Events {
    thrust::device_vector<double> q1, q2, q3, f;
    thrust::device_vector<std::int32_t> n1, n2, n3, cell;
    thrust::device_vector<float> kick;
  } ev_;

  /// Energy the absorptions of a pass gave each partner slot; zero between passes.
  thrust::device_vector<float> dgam_;

  /// Systematic rounding: per-cell running remainders, one array per channel and role
  /// (what becomes of parent 1 (0), of parent 2 (1), the pair count (2)); the scan of
  /// the probabilities; the cell order of unsorted events and their cells alongside.
  struct Rounding {
    std::array<thrust::device_vector<double>, 3 * N_CHANNELS> carry;
    thrust::device_vector<double> cs;
    thrust::device_vector<runko::index_t> perm;
    thrust::device_vector<std::int32_t> keys;
  } round_;

  thrust::device_vector<Staged> staged_;          // products of one channel pass
  thrust::device_vector<std::int32_t> counters_;  // staged count, tail and id counters
  thrust::device_vector<double> dstat_;           // device half of the statistics

  std::uint64_t istep_ { 0 };  // steps taken: the time index of the random streams

  // per-channel slots of dstat_, over the run: one-body events fired, their summed
  // accumulation factor and emitted weight, clamped events, the envelope-weighted rate
  // bound and its weight
  enum DStat : int { d_fired = 0, d_facc, d_wem, d_clamped, d_shnum, d_shden, d_count };

  // slots of counters_: staged records, tail slots per species, products placed in
  // slack, in the tail, and in the tail since an arrival held their slack slot, ids
  // taken per species
  enum Counter : int {
    c_staged = 0,
    c_tail0,
    c_tail1,
    c_tail2,
    c_slack,
    c_tailn,
    c_held,
    c_id0,
    c_id1,
    c_id2,
    c_count
  };

  /// Two-body channel: sides, envelope, sub-steps of strata, thinning and reactions.
  void twobody(Ctx& c, Channel which);

  /// Synchrotron emission of the electrons or the positrons: one lane per slot.
  void synchrotron(Ctx& c, Channel which);

  /// Magnetic pair production of the photons: one lane per slot, then per conversion.
  void breit_wheeler(Ctx& c);

  /// Room for n staged records, none staged yet.
  [[nodiscard]] Stager stager(std::size_t n);

  /// Places the staged products into the slack of their keys or the tail.
  void place(Ctx& c);

  /// Carry array of a role of channel l, seeded on first use and for a new cell count.
  thrust::device_vector<double>& carry(int l, int role, int ncell);

  /// Statistics of channel l: the host half plus its slice of the device read-back h.
  [[nodiscard]] ChannelStats
    channel_snap(std::size_t l, const thrust::host_vector<double>& h) const;

  /// Systematic rounding of the event probabilities q into counts out, on the caller's
  /// queue. sorted: the events come cell-ordered, as the two-body lanes do.
  void select(
    const tyvi::mdgrid_work& w,
    bool sorted,
    int l,
    int role,
    int ncell,
    std::size_t n,
    const std::int32_t* cell,
    const double* q,
    std::int32_t* out);

public:
  Sampler(const Box& b, const Params& p, const std::vector<std::string>& channels);

  [[nodiscard]] std::uint64_t step() const { return istep_; }

  /// Lays out every species (the epoch of the sampler).
  void layout(TileView& t);

  /// Lays out if a species has no valid tables, or on the caller's epoch beat.
  /// The sources and sinks touch the tables, so they call this first.
  void lay_if_needed(TileView& t, bool epoch = false);

  /// One step of every enabled channel, in a random order; lays out where needed.
  void step(TileView& t);

  /// Exact expectation of a two-body channel in one step by a brute-force pair sum:
  /// the events sum_ij P_ij, or the weight they move sum_ij P_ij s_ij.
  [[nodiscard]] double pair_sum(Channel which, const TileView& t, bool weight) const;

  [[nodiscard]] std::map<std::string, double> counts() const;
  [[nodiscard]] std::map<std::string, double> moved() const;
  [[nodiscard]] std::string report() const;

  // ---- sources and sinks of the open box; each accumulates into the budget

  struct [[nodiscard]] Budget {
    double inj_E { 0 }, inj_W { 0 };  // energy and pair depth injected
    double heat_E { 0 };              // energy of the diffusive heating
    double esc_E { 0 }, esc_W { 0 };  // energy and weight escaped
  };
  Budget budget;

  // escape spectrum: NESC log bins of x over [ESC_LO, ESC_HI], NMU linear bins of
  // mu = k_z/x over [-1, 1]
  static constexpr int NESC      = 100;
  static constexpr int NMU       = 32;
  static constexpr double ESC_LO = 1e-8, ESC_HI = 1e2;

  /// Injects pairs at the power ell per step: dN/dgamma ~ gamma^-p on [gmin, gmax],
  /// each lepton isotropic.
  void inject_pairs_plaw(TileView& t, double ell, double p, double gmin, double gmax);

  /// Injects photons of weight w at the power ell per step: a Wien (Gamma(3, Theta))
  /// or Planck spectrum.
  void inject_photons(TileView& t, double ell, double theta, double w, bool planck);

  /// Heats the leptons at the power ell per step by momentum diffusion, D = D_0 p^2.
  void heat_diffusive(TileView& t, double ell);

  /// Removes photons with the slab escape probability at the current pair depth.
  /// tau_T < 0: this tile's own pair depth; dt_esc < 0: the tile is the slab (box.dt).
  void escape_photons(TileView& t, double tau_T = -1.0, double dt_esc = -1.0);

  /// x ell_x of the escaped photons per log bin since the budget was reset.
  [[nodiscard]] std::vector<double> escaped_spectrum() const;

  /// The same per (log x, mu) bin, NESC x NMU row-major.
  [[nodiscard]] std::vector<double> escaped_mu_spectrum() const;

  void reset_budget();

  /// Self-absorption depth across the box of a photon of energy x.
  [[nodiscard]] double tau_ssa(const TileView& t, double x) const;

private:
  thrust::device_vector<double> esc_hist_;     // x ell_x per log bin
  thrust::device_vector<double> esc_mu_hist_;  // the same per (log bin, mu bin)
  thrust::device_vector<double> dbudget_;      // energy and weight of one call
};

}  // namespace qed
