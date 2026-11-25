// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Cell-envelope Monte Carlo of one tile over the energy-binned padded layout of every
// species. Two-body channels (Compton, pair annihilation, photon-photon pair
// production, synchrotron self-absorption) draw per (cell, bin a, bin b) the envelope
// Lambda = c_l N_a N_b 2 sigmahat_ab W_max/N_2 by stratified pairing; one-body channels
// (quantum synchrotron emission, magnetic pair production) run one lane per slot.
// Accumulation merges f = min(w_i/w_j, kappa (gamma-1)/dgamma_max) events of a heavy
// soft photon on a light lepton into one of probability P_ij/f. Products go into the
// slack of their key (reactive at once) or the tail (inert until the next sort).

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

/// Species roles; each container is found from its charge and mass.
enum class Species : std::size_t { ph = 0, em = 1, ep = 2 };

/// What an event does to the partner: kicks it, consumes it, or feeds it energy.
enum class Kind { scatter, convert, absorb };

/// The slot ranges a side spans: the cell's range, one photon bit, or both bits.
enum class Sel { cell, bit0, bit1, both };

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

/// Every channel; the order numbers the random streams of the kernels.
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

// ---------------------------------------------------------------- one-body bookkeeping

/// Accumulation factor of a resolved emission of mean photon energy xm:
/// acc_syn max(1, xm^-2), capped at kappa (gamma - 1)/xm, at least 1.
[[nodiscard]]
constexpr double
  facc_syn(const double gam, double xm, const double acc_syn, const double kappa)
{
  if(acc_syn <= 0.0) { return 1.0; }
  xm = sstd::max(xm, 1e-300);
  return sstd::max(
    sstd::min(acc_syn * sstd::max(1.0, 1.0 / (xm * xm)), kappa * (gam - 1.0) / xm),
    1.0);
}

/// Parent energy after one step of synchrotron cooling above the floor:
/// dgamma/dt = -(2/3)(alpha c/lamC) chi^2 g Phi, in closed form for chi ~ gamma.
/// Phi is the power fraction above the floor, from its rate lam and mean energy xm.
[[nodiscard]]
constexpr double
  gamma_cooled(
    const double gam,
    const double lam,
    const double xm,
    const double b,
    const double N1)
{
  const auto chi  = gam * b;
  const auto g    = syn::g(chi);
  const auto full = (2.0 / 3.0) * chi * chi * g / N1;
  const auto Phi  = full > 0.0 ? lam * xm / sstd::max(full, 1e-300) : 0.0;
  return gam / (1.0 + (2.0 / 3.0) * gam * b * b * g * Phi / N1);
}


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

/// One side of a two-body channel: the one or two sub-ranges of every (cell, bin) seen
/// as one range.
struct [[nodiscard]] SideView {
  using I = runko::index_t;
  const I* base0;     // first slot of sub-range 0, per (cell, bin)
  const I* cnt0;      // its slots (live + appended)
  const I* base1;     // first slot of sub-range 1 (second photon bit), 0 if absent
  const I* cnt1;      // its slots, 0 if absent
  const float* wmax;  // largest weight over both sub-ranges
  const float* emin;  // energy bounds over both sub-ranges
  const float* emax;

  [[nodiscard]] constexpr I N(const int cb) const { return cnt0[cb] + cnt1[cb]; }

  /// Slot of the virtual index v in [0, N) of group cb.
  [[nodiscard]] constexpr I pick(const int cb, const I v) const
  {
    const auto c0 = cnt0[cb];
    return v < c0 ? base0[cb] + v : base1[cb] + v - c0;
  }
};

/// What the sampler needs of the tile for one step.
struct [[nodiscard]] TileView {
  FieldView field;
  std::array<pic::ParticleContainer*, 3> con;  // ph, e-, e+
  std::array<Pop*, 3> pop;                     // their layouts
  std::array<float, 3> mins, maxs;             // tile extent
  std::array<runko::prtc_id_type, 3> id_base;  // first free id per species
  std::array<std::size_t, 3> ids_used;         // returned: ids consumed by the step
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
  /// Decisions taken by systematic rounding: parent updated, kicked or removed; pairs.
  enum class Decision : int { update = 0, kick, remove, pairs, count };

  /// Per-cell running remainders of the rounding, one array per (decision, species).
  std::
    array<thrust::device_vector<double>, 3 * static_cast<std::size_t>(Decision::count)>
      carry_;

  thrust::device_vector<Staged> staged_;          // products of one channel pass
  thrust::device_vector<std::int32_t> counters_;  // staged count, tail and id counters
  thrust::device_vector<double> dstat_;           // device half of the statistics

  /// Energy the absorptions of a pass gave each partner slot; zero between passes.
  thrust::device_vector<float> dgam_;

  // envelope tables per (cell, bin, bin): envelope, rate bound, weight bound; and their
  // sums per (cell, bin)
  thrust::device_vector<double> Lam_, sh_, wm_, rowsum_, colsum_;
  // events per bin pair; the strata bases and per-bin sums of a consumed partner
  thrust::device_vector<std::int32_t> M_, off2_, M2_;
  // first lane of every bin pair, the exclusive scan of M_
  thrust::device_vector<runko::index_t> Mstart_;

  // side tables per (cell, bin) of the incident (1) and partner (2) side: first slot
  // and fill of both sub-ranges, and the bounds over them
  thrust::device_vector<runko::index_t> base1_, cnt1_, base2_, cnt2_;
  thrust::device_vector<float> wmax1_, emin1_, emax1_, wmax2_, emin2_, emax2_;

  // lanes, one per envelope event: the two slots; acc_ lists the accepted lanes
  thrust::device_vector<runko::index_t> lane_i_, lane_j_, acc_;
  thrust::device_vector<std::int32_t> lane_g_;    // bin pair
  thrust::device_vector<double> lane_S_;          // weight the event moves
  thrust::device_vector<std::uint8_t> lane_acc_;  // accepted

  // events, one per accepted lane: decision probabilities (update, kick, lepton count)
  // and their scan; the rounded decisions and the cell
  thrust::device_vector<double> ev_q1_, ev_q2_, ev_q3_, ev_cs_;
  thrust::device_vector<double> ev_f_;  // exchanges the partner takes as one
  thrust::device_vector<std::int32_t> ev_n1_, ev_n2_, ev_n3_, ev_cell_;
  thrust::device_vector<runko::index_t> ev_perm_;  // cell order of unsorted events
  thrust::device_vector<std::int32_t> ev_keys_;    // their cells, sorted alongside
  thrust::device_vector<float> ev_kick_;           // partner kick, 3 per event

  std::uint64_t istep_ { 0 };  // steps taken: the time index of the random streams

  /// Containers, layouts and device views of one step.
  struct Ctx;

  /// Two-body channel l: envelope, sub-cycling, then one pass per sub-step.
  template<typename Ch>
  void twobody(Ctx& c, const Ch& ch, int l);

  /// One pass of a two-body channel over its two sides: count, thin, react, place.
  template<typename Ch>
  void pass(
    Ctx& c,
    const Ch& ch,
    int l,
    int r,
    int nsub,
    const SideView& s1,
    const SideView& s2);

  /// One-body channel: one lane per slot.
  void onebody(Ctx& c, Channel which, int l);

  /// Places the staged products into the slack of their keys or the tail.
  void place(Ctx& c);

  /// Adds the energy absorbed in the na events to the fed leptons, along their momenta.
  void apply_dgam(Ctx& c, Species s, std::size_t na);

  /// Adds the kicks of the na events to the leptons and rewidens the bounds of their
  /// keys.
  void apply_kicks(Ctx& c, Species s, std::size_t na);

  /// Carry array of a decision site, created on first use.
  thrust::device_vector<double>& carry(Decision d, Species s, int ncell);

  /// Statistics of channel l: the host half plus its slice of the device read-back h.
  [[nodiscard]] ChannelStats
    channel_snap(std::size_t l, const thrust::host_vector<double>& h) const;

  /// Systematic rounding of the event probabilities q into counts out, on the caller's
  /// queue. sorted: the events come cell-ordered, as the two-body lanes do.
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
