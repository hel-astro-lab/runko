// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The cell-envelope Monte Carlo of one tile: stratified pairing over the energy-binned padded
// layout of every species, carrying the two-body channels (Compton on either lepton charge,
// pair annihilation, photon-photon pair production, synchrotron self-absorption) and the
// one-body channels (quantum synchrotron emission, magnetic pair production).
//
// The step. Every n_sort steps each species is laid out by its key (reaction cell, energy bin,
// and for photons a random bit) with slack behind every key, and the tables the envelopes read
// are built. The channels then run in a random order, each reading the states the previous one
// left:
//   two-body  per (cell, bin a, bin b) an envelope Lambda = c_l N_a N_b 2 sigmahat_ab W_max/N_2
//             from the tables alone; counts by stochastic rounding, sub-cycled where a bin's
//             envelope exceeds its slot count; two-level stratified pairing on the incident side;
//             the partner side stratified too where an event consumes the partner (conversion),
//             drawn with replacement where it only adds to it (scattering, absorption): the sum of
//             the kicks is what a lepton takes, by an atomic addition. A lane gathers the two
//             weights and tests the weight factor, then the cross section, then the accumulation
//             ratio against a single uniform; the accepted events run the product kernel.
//   one-body  one lane per slot evaluates the rate of its own LP and fires on it.
// Products are staged and then placed in the slack of their own key, reactive at once, or in the
// tail, inert until the next sort; parents are updated in place in the probabilistic
// realization, with per-cell systematic rounding of every decision. The copy number of the
// resolution function is realized exactly in both directions (splitting and Russian roulette).
//
// Accumulation. f = min(w_i/w_j, kappa (gamma-1)/dgamma_max) events of a heavy soft photon on a
// light lepton are merged into one: the pair reacts with P_ij/f, the photon moves the weight f s
// and the lepton takes f exchanges at once; the merged rate is that of a pair of the effective
// weight W_ij = max(w_j, w_i/f_kin), which the maxima of a bin pair bound.

#include "runko/pic/particle.h"
#include "runko/qed/box.h"
#include "runko/qed/field.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/layout.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace qed {

/// The roles the species play; the container of each is found from its charge and mass.
enum class Species : std::size_t { ph = 0, em = 1, ep = 2 };
/// What an event does to the partner: kicks it, consumes it, or feeds it energy.
enum class Kind { scatter, convert, absorb };
/// The slot ranges a side spans: the cell's range, one photon bit, or both bits.
enum class Sel { cell, bit0, bit1, both };

/// The knobs of the sampler.
struct [[nodiscard]] Params {
  int nbin { 16 };          // energy bins per cell and species
  int ncell_edge { 1 };     // reaction sub-cells per tile edge
  int nsort { 10 };         // steps between padded sorts
  int smin { 4 };           // floor on the slack slots of a key
  int n_cap { 32 };         // largest copy number of one reaction
  int nsub_max { 8 };       // sub-steps an additive channel may take before its rows are clamped
  double slack { 0.5 };     // slack slots per key, as a fraction of its live LPs
  double a_res { 0.7 };     // exponent of the photon resolution function
  double x_ref { 1e-5 };    // its anchor: one LP per unit weight there
  double acc_syn { 3.0 };   // coefficient of the emission accumulation, 0 disables
  double kappa { 0.1 };     // merged exchange <= kappa (gamma - 1), 0 disables
  std::uint64_t seed { 1 };
};

/// Every channel, two- and one-body; the order is the kernel numbering of the random streams.
enum class Channel : int {
  compton_m = 0,
  compton_p,
  pairann,
  photann,
  ssa_m,
  ssa_p,
  syn_m,
  syn_p,
  bw,
  count
};
inline constexpr auto N_CHANNELS = static_cast<std::size_t>(Channel::count);

[[nodiscard]] const char* channel_name(Channel c);
[[nodiscard]] bool is_twobody(Channel c);

/// Per-channel counters of a run.
struct [[nodiscard]] ChannelStats {
  double M { 0 };         // envelope lanes (two-body), slots (one-body)
  double virt { 0 };      // two-body lanes on a dead slot
  double wpass { 0 };     // lanes past the weight stage
  double acc { 0 };       // accepted events
  double wmoved { 0 };    // weight moved, sum f s
  double facc { 0 };      // sum of accumulation factors of accepted events
  double wem { 0 };       // emitted weight (one-body emission)
  double clamped { 0 };   // events dropped to keep the strata of a bin distinct
  double shbar_num { 0 }, shbar_den { 0 };  // <sigmahat> of the envelope
  int nsub_max_seen { 1 };  // largest sub-cycle count of a step
  double nsub_steps { 0 };  // steps that needed more than one sub-cycle
};

struct [[nodiscard]] Stats {
  std::size_t steps { 0 };                       // steps taken
  double cell_lanes { 0 };                       // (cell, bin, bin) lanes of the envelope kernels
  std::array<ChannelStats, N_CHANNELS> ch {};    // per channel
  double slackprod { 0 }, tailprod { 0 };              // products placed in slack / in the tail
  std::array<double, 3> slots {}, live {};             // per species: slots, live LPs
};

// ---------------------------------------------------------------- one-body bookkeeping

/// accumulation factor of a resolved emission: coeff max(1, <x>^-2), capped at kappa (gamma-1)/<x>
[[nodiscard]]
constexpr double
  facc_syn(const double gam, double xm, const double acc_syn, const double kappa)
{
  if(acc_syn <= 0.0) { return 1.0; }
  xm = sstd::max(xm, 1e-300);
  return sstd::max(sstd::min(acc_syn * sstd::max(1.0, 1.0 / (xm * xm)), kappa * (gam - 1.0) / xm), 1.0);
}

/// parent energy after one step of continuous synchrotron cooling, radiating above the floor only:
/// dgamma/dt = -(2/3)(alpha c/lamC) chi^2 g Phi integrates in closed form for chi ~ gamma; Phi is the
/// fraction of the power above the floor, from the number rate and the mean energy above it
[[nodiscard]]
constexpr double
  gamma_cooled(const double gam, const double lam, const double xm, const double b, const double N1)
{
  const auto chi  = gam * b;
  const auto g    = syn::g(chi);
  const auto full = (2.0 / 3.0) * chi * chi * g / N1;
  const auto Phi  = full > 0.0 ? lam * xm / sstd::max(full, 1e-300) : 0.0;
  return gam / (1.0 + (2.0 / 3.0) * gam * b * b * g * Phi / N1);
}


/// A staged product: momentum, position, weight per copy, copy count, destination, species.
struct [[nodiscard]] Staged {
  float u[3];            // momentum
  float pos[3];          // position (the parent's)
  float w;               // weight of each copy
  std::int32_t n_r;      // copies to place
  std::int32_t cell;     // reaction cell
  std::int32_t species;  // 0 ph, 1 e-, 2 e+
  std::int32_t bit;      // photon sub-range bit
};

/// One side of a two-body channel: the one or two sub-ranges of every (cell, bin) as one virtual range.
struct [[nodiscard]] SideView {
  using I = runko::index_t;
  const I* base0;      // first slot of sub-range 0 of every (cell, bin)
  const I* cnt0;       // its slot count (live + appended)
  const I* base1;      // first slot of sub-range 1 (the second photon bit), 0 if absent
  const I* cnt1;       // its count, 0 if absent
  const float* wmax;   // largest weight over both sub-ranges
  const float* emin;   // energy bounds over both sub-ranges
  const float* emax;
  [[nodiscard]] constexpr I N(const int cb) const { return cnt0[cb] + cnt1[cb]; }
  /// virtual index v in [0, N) of group cb -> slot
  [[nodiscard]] constexpr I pick(const int cb, const I v) const
  {
    const auto c0 = cnt0[cb];
    return v < c0 ? base0[cb] + v : base1[cb] + v - c0;
  }
};

/// What the sampler needs of the tile for one step.
struct [[nodiscard]] TileView {
  FieldView field;                              // the field the channels see
  std::array<pic::ParticleContainer*, 3> con;  // ph, e-, e+
  std::array<Pop*, 3> pop;                     // their layouts
  std::array<float, 3> mins, maxs;             // tile extent
  std::array<runko::prtc_id_type, 3> id_base;   // first free id of every species
  std::array<std::size_t, 3> ids_used;          // returned: ids consumed by the step
};

/// The sampler state of a tile: scratch buffers, carries and counters that persist over steps.
class [[nodiscard]] Sampler {
public:
  Box box;                                   // physical normalization
  Params par;                                // the knobs
  std::array<bool, N_CHANNELS> enabled {};   // channels to run
  Stats stats;                               // counters of the run

private:
  /// the decisions taken by systematic rounding: a parent updated, kicked or removed, pairs made
  enum class Decision : int { update = 0, kick, remove, pairs, count };
  // per-cell running remainders of the rounding, one array per (decision, species) site
  std::array<thrust::device_vector<double>, 3 * static_cast<std::size_t>(Decision::count)> carry_;
  // staged products of one channel pass
  thrust::device_vector<Staged> staged_;
  thrust::device_vector<std::int32_t> counters_;  // staged count, tail counters, id counters
  thrust::device_vector<double> dstat_;           // per-channel statistics accumulated on the device over the run
  // the energy the absorptions of a pass gave each partner slot; zero between passes
  thrust::device_vector<float> dgam_;
  // envelope tables of a pass, per (cell, bin, bin) pair unless noted
  thrust::device_vector<double> Lam_, sh_, wm_, rowsum_, colsum_;  // envelope, rate bound, weight bound; sums per (cell, bin)
  thrust::device_vector<std::int32_t> M_, off1_, off2_, M1_, M2_, sh1_, sh2_, kk_;  // events, strata bases, per-bin sums and shifts, stratum offsets
  thrust::device_vector<runko::index_t> Mstart_;  // first lane of every pair (exclusive scan of M_)
  // side tables, per (cell, bin) of the incident (1) and the partner (2) side
  thrust::device_vector<runko::index_t> base1_, cnt1_, base2_, cnt2_;  // first slot and fill of the two sub-ranges
  thrust::device_vector<float> wmax1_, emin1_, emax1_, wmax2_, emin2_, emax2_;  // bounds over the sub-ranges
  // lanes (one per envelope event) and events (one per accepted lane)
  thrust::device_vector<runko::index_t> lane_i_, lane_j_, acc_;  // the two slots of a lane; the accepted lanes
  thrust::device_vector<std::int32_t> lane_g_;                   // the bin pair of a lane
  thrust::device_vector<double> lane_f_;                         // accumulation factor of a lane
  thrust::device_vector<std::uint8_t> lane_acc_;                 // lane accepted
  thrust::device_vector<double> ev_q1_, ev_q2_, ev_q3_, ev_cs_;  // decision probabilities (update/kick/lepton count) and their scan
  thrust::device_vector<double> ev_fs_;                          // weight moved per event, f s
  thrust::device_vector<std::int32_t> ev_n1_, ev_n2_, ev_n3_, ev_cell_;  // the rounded decisions and the cell
  thrust::device_vector<runko::index_t> ev_perm_;                // cell-sorted order of unsorted events
  thrust::device_vector<std::int32_t> ev_keys_;                  // their cells, sorted alongside
  thrust::device_vector<float> ev_p_;  // product momenta, 6 per event
  std::uint64_t istep_ { 0 };          // steps taken: the time index of the random streams

  struct Ctx;  // the containers, layouts and device views of one step
  /// a two-body channel: envelope, sub-cycling, then one pass per sub-step
  template<typename Ch>
  void twobody(Ctx& c, const Ch& ch, int l);
  /// one pass of a two-body channel over its two sides: count, thin, react, place
  template<typename Ch>
  void pass(Ctx& c, const Ch& ch, int l, int r, int nsub, const SideView& s1, const SideView& s2);
  /// a one-body channel: one lane per slot
  void onebody(Ctx& c, Channel which, int l);
  /// the staged products into the slack of their keys or the tail
  void place(Ctx& c);
  /// the summed absorbed energies onto the fed leptons, along their momenta, over the na events
  void apply_dgam(Ctx& c, Species s, std::size_t na);
  /// rewiden the bounds of the keys of the kicked leptons, over the na events
  void touch_kicked(Ctx& c, Species s, std::size_t na);
  /// the carry array of a decision site, created on first use
  thrust::device_vector<double>& carry(Decision d, Species s, int ncell);
  /// a channel's statistics: the host half plus the device slice of a read-back copy
  [[nodiscard]] ChannelStats channel_snap(std::size_t l, const thrust::host_vector<double>& h) const;
  /// systematic rounding of the probabilities q of the events into counts out, on the
  /// caller's queue (stream-ordered with the caller's kernels); sorted: the events come
  /// cell-ordered, which the two-body lanes are by construction
  void select(
    const tyvi::mdgrid_work& w,
    bool sorted,
    Decision d,
    Species s,
    int ncell,
    std::size_t n,
    const std::int32_t* cell,
    const double* q,
    std::int32_t* out);

public:
  Sampler(const Box& b, const Params& p, const std::vector<std::string>& channels);

  [[nodiscard]] std::uint64_t step() const { return istep_; }

  /// The layout of every species (the epoch of the sampler).
  void layout(TileView& t);

  /// Lay out when a species has no valid tables (or on the caller's epoch beat): the sources
  /// and sinks touch the tables, so they call this first.
  void lay_if_needed(TileView& t, bool epoch = false);

  /// One step of every enabled channel, in a random order; lays out where needed.
  void step(TileView& t);

  /// Exact expectation of a two-body channel in one step by a brute-force pair sum: the events
  /// sum_ij P_ij, or the weight they move sum_ij P_ij s_ij.
  [[nodiscard]] double pair_sum(Channel which, const TileView& t, bool weight) const;

  [[nodiscard]] std::map<std::string, double> counts() const;
  [[nodiscard]] std::map<std::string, double> moved() const;
  [[nodiscard]] std::string report() const;

  // ---- sources and sinks of the open box; every one accumulates into the budget
  struct [[nodiscard]] Budget {
    double inj_E { 0 }, inj_W { 0 };   // energy and pair depth injected
    double heat_E { 0 };               // energy of the diffusive heating
    double esc_E { 0 }, esc_W { 0 };   // energy and weight escaped
  };
  Budget budget;
  static constexpr int NESC = 100;     // log bins of the escaping spectrum over [1e-8, 1e2]
  static constexpr int NMU  = 32;      // linear bins of the escape direction mu = k_z/x over [-1, 1]
  static constexpr double ESC_LO = 1e-8, ESC_HI = 1e2;

  /// pairs at the power ell per step: dN/dgamma ~ gamma^-p on [gmin, gmax], each lepton isotropic
  void inject_pairs_plaw(TileView& t, double ell, double p, double gmin, double gmax);
  /// photons at the power ell per step, weight w each: a Wien (Gamma(3, Theta)) or Planck spectrum
  void inject_photons(TileView& t, double ell, double theta, double w, bool planck);
  /// stochastic heating of the leptons at the power ell per step: momentum diffusion with D = D_0 p^2
  void heat_diffusive(TileView& t, double ell);
  /// remove photons with the slab escape probability at the current pair depth
  /// tau_T < 0: this tile's own pair depth; dt_esc < 0: the tile is the slab (box.dt)
  void escape_photons(TileView& t, double tau_T = -1.0, double dt_esc = -1.0);
  /// x ell_x of the escaped photons per log bin, summed since the budget was reset
  [[nodiscard]] std::vector<double> escaped_spectrum() const;
  /// x ell_x of the escaped photons per (log x, mu) bin, NESC x NMU row-major, since the reset
  [[nodiscard]] std::vector<double> escaped_mu_spectrum() const;
  void reset_budget();
  /// self-absorption depth across the box of a photon of energy x against the current leptons
  [[nodiscard]] double tau_ssa(const TileView& t, double x) const;

private:
  thrust::device_vector<double> esc_hist_;      // NESC: x ell_x per log bin
  thrust::device_vector<double> esc_mu_hist_;   // NESC x NMU: the same per escape direction
  thrust::device_vector<double> dbudget_;   // energy and weight accumulated by one source or sink
};

}  // namespace qed
