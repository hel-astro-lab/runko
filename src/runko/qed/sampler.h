// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

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
struct Params {
  /// Energy bins per cell and species.
  int nbin { 16 };
  /// Reaction sub-cells per tile edge.
  int ncell_edge { 1 };
  /// Steps between padded sorts.
  int nsort { 10 };
  /// Floor on the slack slots of a key.
  int smin { 4 };
  /// Largest copy number of one reaction.
  int n_cap { 32 };
  /// Sub-steps of an additive channel before rows are clamped.
  int nsub_max { 8 };
  /// Slack slots per key, as a fraction of its live LPs.
  double slack { 0.5 };
  /// Exponent of the photon resolution function.
  double a_res { 0.7 };
  /// Anchor of the photon resolution function: one LP per unit weight there.
  double x_ref { 1e-5 };
  /// Coefficient of the emission accumulation; 0 disables it.
  double acc_syn { 3.0 };
  /// Bound kappa (gamma - 1) on a merged exchange; 0 disables merging.
  double kappa { 0.1 };
  std::uint64_t seed { 1 };
};

/// Per-channel counters of a run.
struct ChannelStats {
  /// Envelope lanes (two-body) or slots (one-body).
  double M { 0 };
  /// Two-body lanes on a dead slot.
  double virt { 0 };
  /// Lanes past the weight stage.
  double wpass { 0 };
  /// Accepted events.
  double acc { 0 };
  /// Weight moved: the summed sizes of the accepted events.
  double wmoved { 0 };
  /// Summed accumulation factors of the accepted events.
  double facc { 0 };
  /// Emitted weight (one-body emission).
  double wem { 0 };
  /// Events dropped to keep the strata of a bin distinct.
  double clamped { 0 };
  /// Numerator and denominator of the mean sigmahat of the envelope.
  double shbar_num { 0 }, shbar_den { 0 };
  /// Largest sub-cycle count of a step.
  int nsub_max_seen { 1 };
  /// Steps with more than one sub-cycle.
  double nsub_steps { 0 };
};

/// Counters of a run.
struct Stats {
  std::size_t steps { 0 };
  /// Lanes (cell, bin, bin) of the envelope kernels.
  double cell_lanes { 0 };
  std::array<ChannelStats, n_channels> ch {};
  /// Products placed in the slack of their key and in the tail.
  double slackprod { 0 }, tailprod { 0 };
  /// Tail products whose slack slot held an arrival.
  double heldprod { 0 };
  /// Slots and live particles per species.
  std::array<double, 3> slots {}, live {};
};

/// A product of a channel pass, staged before placement.
struct Staged {
  float u[3];
  /// The parent's position.
  float pos[3];
  /// Weight of each copy.
  float w;
  /// Copies to place.
  std::int32_t n_r;
  std::int32_t cell;
  /// 0 photon, 1 electron, 2 positron.
  std::int32_t species;
  /// Photon sub-range bit.
  std::int32_t bit;
};

/// Appends product records to the staging buffer from inside a kernel.
struct Stager {
  Staged* rec;
  /// Records staged so far.
  std::int32_t* count;

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

/// Cell-envelope Monte Carlo of one tile over the energy-binned padded layout.
///
/// Two-body channels (Compton, pair annihilation, photon-photon pair production,
/// synchrotron self-absorption) draw per (cell, bin a, bin b) the envelope
/// Lambda = c_l N_a N_b 2 sigmahat_ab W_max/N_2 by stratified pairing (twobody.c++);
/// one-body channels (quantum synchrotron emission, magnetic pair production) run one
/// lane per slot (onebody.c++). What an event moves and how often it fires is the rule
/// of weights.h. Products go into the slack of their key (reactive at once) or the
/// tail (inert until the next layout). The scratch buffers are members to keep their
/// allocations between passes; the rounding carries and the counters persist over
/// steps.
class [[nodiscard]] Sampler {
public:
  Box box;
  Params par;
  /// Indexed by Channel.
  std::array<bool, n_channels> enabled {};
  Stats stats;

  /// What the sources and sinks of the open box moved since the last reset.
  struct Budget {
    /// Energy and pair depth injected.
    double inj_E { 0 }, inj_W { 0 };
    /// Energy of the diffusive heating.
    double heat_E { 0 };
    /// Energy and weight escaped.
    double esc_E { 0 }, esc_W { 0 };
  };
  Budget budget;

  /// Escape spectrum: n_esc log bins of x over [esc_lo, esc_hi], n_mu linear bins of
  /// mu = k_z/x over [-1, 1].
  static constexpr int n_esc      = 100;
  static constexpr int n_mu       = 32;
  static constexpr double esc_lo = 1e-8, esc_hi = 1e2;

private:
  template<typename Ch>
  friend class TwoBody;

  /// One side of a two-body channel, as of the start of the channel.
  ///
  /// Per (cell, bin): the first slot, the fill of its one or two sub-ranges (the
  /// second at ncb + cb) and the bounds over them.
  struct SideTable {
    thrust::device_vector<runko::index_t> base, cnt;
    thrust::device_vector<float> wmax, emin, emax;
  };
  /// Incident and partner side.
  std::array<SideTable, 2> sides_;

  /// Envelope tables per (cell, bin a, bin b).
  ///
  /// The envelope, the rate bound and the weight bound; the envelope summed over b
  /// (row) and over a (col).
  struct Envelope {
    thrust::device_vector<double> Lam, bound, wmax, row, col;
  } env_;

  /// Envelope events of a sub-step per (cell, bin a, bin b).
  ///
  /// The count M and its exclusive scan, the first lane; for a consumed partner the
  /// strata base along a and the events per partner bin.
  struct Strata {
    thrust::device_vector<std::int32_t> M, off2, M2;
    thrust::device_vector<runko::index_t> first;
  } strata_;

  /// How the thinning of a lane ended.
  ///
  /// A lane on an unkeyed slot counts as dead; weight_passed is a lane past the
  /// weight stage and rejected by the rate.
  enum class Outcome : std::uint8_t { weight_failed, accepted, dead, weight_passed };

  /// One envelope event.
  ///
  /// The slot pair i, j it proposes, the weight S its event would move, its bin pair
  /// g and how its thinning ended.
  struct Lane {
    runko::index_t i, j;
    double S;
    std::int32_t g;
    Outcome outcome;
  };
  thrust::device_vector<Lane> lanes_;

  /// Accepted lanes (two-body) or fired slots (one-body).
  thrust::device_vector<runko::index_t> acc_;
  /// One-body: whether the slot fired.
  thrust::device_vector<std::uint8_t> fired_;

  /// Accepted events.
  ///
  /// The probabilities q1, q2 that parent 1, 2 is updated or removed and the pair
  /// count q3, their rounded decisions n1, n2, n3, the cell, the exchanges f the
  /// partner takes as one and its momentum kick (3 per event).
  struct Events {
    thrust::device_vector<double> q1, q2, q3, f;
    thrust::device_vector<std::int32_t> n1, n2, n3, cell;
    thrust::device_vector<float> kick;
  } ev_;

  /// Energy the absorptions of a pass gave each partner slot; zero between passes.
  thrust::device_vector<float> dgam_;

  /// State of the systematic rounding.
  ///
  /// Per-cell running remainders, one array per channel and role (parent 1 (0),
  /// parent 2 (1), the pair count (2)); the scan of the probabilities; the cell order
  /// of unsorted events and their cells alongside.
  struct Rounding {
    std::array<thrust::device_vector<double>, 3 * n_channels> carry;
    thrust::device_vector<double> cs;
    thrust::device_vector<runko::index_t> perm;
    thrust::device_vector<std::int32_t> keys;
  } round_;

  /// Products of one channel pass.
  thrust::device_vector<Staged> staged_;
  /// Their count, one element.
  thrust::device_vector<std::int32_t> nstaged_;

  /// Counters of one placement, summed on the device in the one element of placed_.
  struct Placed {
    /// Tail slots taken, per species.
    std::int32_t tail[3];
    /// Ids taken, per species.
    std::int32_t id[3];
    /// Products placed in the slack of their key.
    std::int32_t slack;
    /// Products placed in the tail.
    std::int32_t tailn;
    /// Tail products whose slack slot held an arrival.
    std::int32_t held;
  };
  thrust::device_vector<Placed> placed_;

  /// Per-channel statistics the kernels sum on the device over the run.
  ///
  /// The one-body events fired, their summed accumulation factor and emitted weight;
  /// the clamped two-body events, the envelope-weighted rate bound and its weight.
  struct DeviceStats {
    double fired, facc, wem, clamped, shnum, shden;
  };
  thrust::device_vector<DeviceStats> dstat_;

  /// Escaped x ell_x per log bin.
  thrust::device_vector<double> esc_hist_;
  /// Escaped x ell_x per (log bin, mu bin).
  thrust::device_vector<double> esc_mu_hist_;

  /// Energy and weight a source or sink moved in one call, summed on the device in
  /// the one element of tally_.
  struct Tally {
    double E, W;
  };
  thrust::device_vector<Tally> tally_;

  /// Steps taken: the time index of the random streams.
  std::uint64_t istep_ { 0 };

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

  /// Statistics of channel l: the host half plus its part of the device read-back h.
  [[nodiscard]] ChannelStats
    channel_snap(std::size_t l, const thrust::host_vector<DeviceStats>& h) const;

  /// Systematic rounding of the event probabilities q into counts out, on the
  /// caller's queue.
  ///
  /// Set sorted when the events come cell-ordered, as the two-body lanes do.
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
  ///
  /// The sources and sinks touch the tables, so they call this first.
  void lay_if_needed(TileView& t, bool epoch = false);

  /// One step of every enabled channel, in a random order; lays out where needed.
  void step(TileView& t);

  /// Exact expectation of a two-body channel in one step by a brute-force pair sum.
  ///
  /// The events sum_ij P_ij, or with weight set the weight they move,
  /// sum_ij P_ij s_ij.
  [[nodiscard]] double pair_sum(Channel which, const TileView& t, bool weight) const;

  [[nodiscard]] std::map<std::string, double> counts() const;
  [[nodiscard]] std::map<std::string, double> moved() const;
  [[nodiscard]] std::string report() const;

  /// Injects pairs at the power ell per step: dN/dgamma ~ gamma^-p on [gmin, gmax],
  /// each lepton isotropic.
  void inject_pairs_plaw(TileView& t, double ell, double p, double gmin, double gmax);

  /// Injects photons of weight w at the power ell per step: a Wien (Gamma(3, Theta))
  /// or Planck spectrum.
  void inject_photons(TileView& t, double ell, double theta, double w, bool planck);

  /// Heats the leptons at the power ell per step by momentum diffusion, D = D_0 p^2.
  void heat_diffusive(TileView& t, double ell);

  /// Removes photons with the slab escape probability at the current pair depth.
  ///
  /// tau_T < 0 takes this tile's own pair depth; dt_esc < 0 takes the tile as the
  /// slab (box.dt).
  void escape_photons(TileView& t, double tau_T = -1.0, double dt_esc = -1.0);

  /// Escaped x ell_x per log bin since the budget was reset.
  [[nodiscard]] std::vector<double> escaped_spectrum() const;

  /// The same per (log x, mu) bin, n_esc x n_mu row-major.
  [[nodiscard]] std::vector<double> escaped_mu_spectrum() const;

  /// Starts a new accounting window: the budget and the escape spectra.
  void reset_budget();

  /// Self-absorption depth across the box of a photon of energy x.
  [[nodiscard]] double tau_ssa(const TileView& t, double x) const;
};

}  // namespace qed
