// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/sampler.h"

#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "thrust/binary_search.h"
#include "thrust/copy.h"
#include "thrust/fill.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/reduce.h"
#include "thrust/scan.h"
#include "thrust/sort.h"
#include "thrust/transform.h"
#include "thrust/transform_reduce.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <format>
#include <sstream>
#include <stdexcept>

namespace qed {

using I = runko::index_t;
namespace {

// channel names in the order of the Channel enum, the keys of counts() and moved()
constexpr const char* NAMES[] = { "compton-", "compton+", "pairann", "photann", "ssa-",
                                  "ssa+",     "syn-",     "syn+",    "bw" };

// ---------------------------------------------------------------- channel records
// The physics bound to records, nothing dispatched on names later: the sampler kernels are
// templates over these types.
// Compton: photon on a lepton of either charge; scatters, bounded by sigma_KN at the top of the bins
struct ChCompton {
  Species t2;                                          // the lepton species, e- or e+
  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::scatter;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  static constexpr double g_min = 1.0;                 // every lepton scatters
  constexpr double rate(const Vec3& k, const Vec3& u) const { return comp::rate(k, u); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  {
    return comp::bound(lo1, hi1, lo2, hi2);
  }
  constexpr double dgmax(double x, double g) const { return comp::dgmax(x, g); }
};

// pair annihilation: e- on e+; converts both into photons, bounded by the annihilation bracket
struct ChAnn {
  static constexpr Species t1   = Species::em;
  static constexpr Species t2   = Species::ep;
  static constexpr Kind kind    = Kind::convert;
  static constexpr Sel sel1     = Sel::cell;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = false;
  constexpr double rate(const Vec3& um, const Vec3& up) const { return ann::rate(um, up); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  {
    return ann::bound(lo1, hi1, lo2, hi2);
  }
};

// photon-photon pair production: photon bit 0 on photon bit 1; converts, bounded by sigma_gg
struct ChGG {
  static constexpr Species t1   = Species::ph;
  static constexpr Species t2   = Species::ph;
  static constexpr Kind kind    = Kind::convert;
  static constexpr Sel sel1     = Sel::bit0;
  static constexpr Sel sel2     = Sel::bit1;
  static constexpr double c_l   = 2.0;  // the bit-split same-species channel
  static constexpr bool has_dgm = false;
  constexpr double rate(const Vec3& k1, const Vec3& k2) const { return gg::rate(k1, k2); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  {
    return gg::bound(lo1, hi1, lo2, hi2);
  }
};

// synchrotron self-absorption: photon on a lepton of either charge; the lepton absorbs it,
// bounded by a scan of the rate over the lepton bin at the photon bin's floor
struct ChSSA {
  Species t2;                                          // the lepton species, e- or e+
  double b_tile;                                       // the tile field, an upper bound of every lepton's b_ssa
  double x_min;                                        // floor of the synchrotron band
  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::absorb;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  static constexpr double g_min = ssa::GMIN;           // no absorption below the threshold
  /// the rate at the absorbing lepton's field b (its pitch angle in the tile field)
  constexpr double rate(const Vec3& k, const Vec3& u, const double b) const
  {
    return ssa::rate(k, u, b, x_min);
  }
  /// the bound at the tile field, which every lepton's b_ssa is below
  constexpr double bound(double lo1, double /*hi1*/, double lo2, double hi2) const
  {
    return ssa::bound(b_tile, x_min, lo1, lo2, hi2);
  }
  constexpr double dgmax(double x, double) const { return ssa::dgmax(x, 0.0); }
};

// ---------------------------------------------------------------- device views
// the mdspan views a container hands out, captured by value into the kernels
using VecMds    = decltype(std::declval<pic::ParticleContainer&>().vel_mds());
using IdMds     = decltype(std::declval<pic::ParticleContainer&>().ids_mds());
using WgtMds    = decltype(std::declval<pic::ParticleContainer&>().wgt_mds());

/// one species as the kernels see it
struct PView {
  VecMds pos, vel;   // positions (cells) and momenta (m_e c; k for photons)
  IdMds ids;         // dead_prtc_id marks a free slot
  WgtMds wgt;        // per-particle weights, empty for leptons
  bool weighted;     // whether wgt is used (photons) or w0 (leptons)
  float w0;          // the common weight of an unweighted species
  float mass;        // 0 for photons, 1 for leptons
  LayoutView lay;    // the key tables of the species, for touch()

  [[nodiscard]] constexpr bool alive(const I i) const { return ids[i][] != runko::dead_prtc_id; }
  // momentum, weight and energy of slot i as doubles
  [[nodiscard]] constexpr Vec3 u(const I i) const
  {
    return Vec3(static_cast<double>(vel[i][0]), static_cast<double>(vel[i][1]),
                static_cast<double>(vel[i][2]));
  }
  [[nodiscard]] constexpr double w(const I i) const
  {
    return static_cast<double>(weighted ? wgt[i][] : w0);
  }
  [[nodiscard]] constexpr double e(const I i) const
  {
    const auto v = u(i);
    return sstd::sqrt(static_cast<double>(mass) * mass + toolbox::dot(v, v));
  }
  // write back a momentum or a weight; kill frees the slot
  constexpr void set_u(const I i, const Vec3& v) const
  {
    for(auto k = 0uz; k < 3uz; ++k) { vel[i][k] = static_cast<float>(v[k]); }
  }
  constexpr void set_w(const I i, const double v) const
  {
    if(weighted) { wgt[i][] = static_cast<float>(v); }
  }
  constexpr void kill(const I i) const
  {
    ids[i][] = runko::dead_prtc_id;
    if(weighted) { wgt[i][] = 0.0f; }
  }
  /// a slot updated in place: widen the bounds of its key, or, if the kick took it out of its
  /// energy bin, drop it from the key -- inert until the next layout, like a tail product --
  /// so that the bounds of a bin stay those of the bin
  constexpr void touch(const I i) const
  {
    const auto k = lay.skey[i];
    if(k < 0) { return; }
    const auto ee = static_cast<float>(e(i));
    if(lay.bin_of(ee) != lay.bin_of_key(k)) {
      lay.skey[i] = -1;
      return;
    }
    const auto ww = static_cast<float>(w(i));
    sstd::atomic_max(&lay.emax[k], ee);
    sstd::atomic_min(&lay.emin[k], ee);
    sstd::atomic_max(&lay.wmax[k], ww);
  }
};

// the view of a container alone: what the sources, the sinks and the pair sum read
PView
  make_view(pic::ParticleContainer& con)
{
  return PView { .pos      = con.pos_mds(),
                 .vel      = con.vel_mds(),
                 .ids      = con.ids_mds(),
                 .wgt      = con.wgt_mds(),
                 .weighted = con.has_particle_weights(),
                 .w0       = static_cast<float>(con.weight()),
                 .mass     = static_cast<float>(con.args().mass),
                 .lay      = {} };
}

// the view of a container and its layout, refreshed after every resize
PView
  make_view(pic::ParticleContainer& con, Pop& pop)
{
  pop.follow_container(con);
  auto v = make_view(con);
  v.lay  = pop.view();
  return v;
}

/// the rate of the pair (p1, p2); absorption takes the partner's field across its pitch angle
template<typename Ch>
[[nodiscard]]
constexpr double
  pair_rate(const Ch& ch, const FieldView& fv, const Vec3& p1, const Vec3& p2)
{
  if constexpr(Ch::kind == Kind::absorb) {
    return ch.rate(p1, p2, fv.b_ssa(p2 / sstd::max(toolbox::norm(p2), 1e-300)));
  } else {
    return ch.rate(p1, p2);
  }
}

/// offset inside a range of N slots cut into M strata of a slot drawn uniformly from stratum p
[[nodiscard]]
constexpr I
  stratum(const std::uint64_t p, const std::uint64_t M, const std::uint64_t N, runko::Rng& rng)
{
  const auto B0 = (p * N) / M;        // first slot of stratum p
  const auto B1 = ((p + 1) * N) / M;  // one past its last
  return static_cast<I>(B0 + static_cast<std::uint64_t>(rng.uniform() * static_cast<double>(B1 - B0)));
}

/// index of the last start <= v in a sorted array of n starts
[[nodiscard]]
constexpr int
  segment_of(const I* start, const int n, const I v)
{
  auto lo = 0, hi = n;  // start[lo] <= v < start[hi]
  while(hi - lo > 1) {
    const auto mid = (lo + hi) / 2;
    if(start[mid] <= v) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

/// resolution function of the photons: LPs per unit weight go as (x/x_ref)^a
[[nodiscard]]
constexpr double
  f_res(const double e, const double a, const double x_ref)
{
  return sstd::pow(e / x_ref, a);
}

/// realize the copy number n of a product of weight s: (n_r, weight of each copy). n >= 1 splits
/// into floor(n) or floor(n) + 1 copies of s/n_r; n < 1 keeps one copy with probability n at s/n
/// (Russian roulette). Both keep <n_r> = min(n, n_cap).
struct Copies {
  int n_r;
  double w_each;
};
[[nodiscard]]
constexpr Copies
  copies(double n, const double s, runko::Rng& rng, const int n_cap)
{
  n              = sstd::min(n, static_cast<double>(n_cap));
  const auto n_r = static_cast<int>(sstd::floor(n + rng.uniform()));
  const auto den = n >= 1.0 ? static_cast<double>(sstd::max(n_r, 1)) : sstd::max(n, 1e-300);
  return Copies { n_r, s / den };
}

/// accumulation factor: events a lepton can take at once before its state moves by kappa (gamma-1)
template<typename Ch>
[[nodiscard]]
constexpr double
  fkin(const Ch& ch, const double kappa, const double x, const double g)
{
  if constexpr(not Ch::has_dgm) {
    return 1.0;
  } else {
    if(kappa <= 0.0) { return 1.0; }
    return sstd::max(1.0, kappa * (g - 1.0) / sstd::max(ch.dgmax(x, g), 1e-300));
  }
}

/// a staged record of one product at a given position and cell
[[nodiscard]]
constexpr Staged
  staged_at(const Vec3& u, const float* pos, const double w, const int n_r, const int cell, const int species, const int bit)
{
  return Staged { .u       = { static_cast<float>(u[0]), static_cast<float>(u[1]), static_cast<float>(u[2]) },
                  .pos     = { pos[0], pos[1], pos[2] },
                  .w       = static_cast<float>(w),
                  .n_r     = n_r,
                  .cell    = cell,
                  .species = species,
                  .bit     = bit };
}

// slots of the per-pass statistics accumulated on the device (dstat_): lanes on a dead slot,
// lanes past the weight stage, moved weight, summed accumulation factor, rejection proposals
// and their count, clamped events, emitted weight
enum DStat : int { d_virt = 0, d_wpass, d_wmoved, d_facc, d_trials, d_trials_n, d_clamped, d_wem, d_count };
// slots of the device counters (counters_): staged records, tail slots taken per species,
// products placed in slack and in the tail, ids consumed per species
enum Counter : int { c_staged = 0, c_tail0, c_tail1, c_tail2, c_slack, c_tailn, c_id0, c_id1, c_id2, c_count };

}  // namespace

// ---------------------------------------------------------------- the per-step context
struct Sampler::Ctx {
  TileView* t;                        // the tile's containers, layouts, extent and ids
  std::array<PView, 3> v;             // views of ph, e-, e+; stale after a resize
  int ncell;                          // reaction cells of the tile
  int nb;                             // energy bins per cell

  Ctx(TileView& tile, const int nbin) : t { &tile }, v {}, ncell { 0 }, nb { nbin }
  {
    refresh_all();
    ncell = v[0].lay.ncell;
  }

  // rebuild the views after a container was resized
  void refresh_all()
  {
    for(auto i = 0uz; i < 3uz; ++i) { v[i] = make_view(*t->con[i], *t->pop[i]); }
  }
};

const char*
  channel_name(const Channel c)
{
  return NAMES[static_cast<int>(c)];
}

bool
  is_twobody(const Channel c)
{
  return static_cast<int>(c) <= static_cast<int>(Channel::ssa_p);
}

Sampler::Sampler(const Box& b, const Params& p, const std::vector<std::string>& channels) :
  box { b }, par { p }
{
  // enable the channels: the field channels need b_perp > 0, a given list restricts further
  const auto field = b.b_perp() > 0.0;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    const auto c    = static_cast<Channel>(i);
    const auto need = static_cast<int>(c) >= static_cast<int>(Channel::ssa_m);  // needs a field
    auto on         = not need or field;
    if(not channels.empty()) {
      on = on and std::ranges::find(channels, std::string(NAMES[i])) != channels.end();
    }
    enabled[i] = on;
  }
  // the small device accumulators every pass zeroes and reads back
  counters_.assign(c_count, 0);
  dstat_.assign(d_count, 0.0);
  dbudget_.assign(2, 0.0);
}

// the per-cell remainders of one decision site, created on first use
thrust::device_vector<double>&
  Sampler::carry(const std::string& site, const int ncell)
{
  auto it = carry_.find(site);
  if(it == carry_.end()) {
    // seeded from the stream of the site, so that a run is reproducible
    auto h   = thrust::host_vector<double>(static_cast<std::size_t>(ncell));
    auto key = std::uint64_t { 0 };
    for(const auto ch: site) { key = runko::mix64(key ^ static_cast<std::uint64_t>(ch)); }
    auto rng = runko::Rng { runko::rng_key(par.seed, 0, 0, key), 0 };
    for(auto& x: h) { x = rng.uniform(); }
    it = carry_.emplace(site, thrust::device_vector<double>(h)).first;
  }
  return it->second;
}

// Every per-event decision -- a copy count, an update, a removal -- is taken against a running
// remainder of its cell rather than an independent uniform: the marginal probability of each
// event is unchanged, and the summed error of a cell is bounded by one instead of growing as a
// random walk. Events come cell-sorted, so the running sum is a segmented scan.
void
  Sampler::select(
    const tyvi::mdgrid_work& w,
    const bool sorted,
    const std::string& site,
    const int ncell,
    const std::size_t n,
    const std::int32_t* cell,
    const double* q,
    std::int32_t* out)
{
  if(n == 0) { return; }
  auto& car        = carry(site, ncell);
  const auto car_p = thrust::raw_pointer_cast(car.data());   // remainder per cell
  ev_cs_.resize(n);
  const auto cs = thrust::raw_pointer_cast(ev_cs_.data());     // running sum of q within the cell
  const auto ord0 = thrust::counting_iterator<I>(0);
  // the base of a cell is read by every event of the cell and updated by all of them, so the
  // update is a second pass: reading and adding in one kernel would race
  if(sorted) {
    // running sum of the probabilities within each cell
    thrust::inclusive_scan_by_key(
      w.on_this(), cell, cell + static_cast<std::ptrdiff_t>(n), q, cs);
    // one lane per event: its count from the remainder and the running sum
    thrust::for_each(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      const auto base = car_p[cell[i]];
      const auto m    = sstd::floor(base + cs[i]) - sstd::floor(base + cs[i] - q[i]);
      out[i]          = static_cast<std::int32_t>(m);
    });
  } else {
    // the tail of a species is not cell-ordered: sort a permutation by cell first
    ev_perm_.resize(n);
    ev_keys_.resize(n);
    thrust::copy(w.on_this(), cell, cell + static_cast<std::ptrdiff_t>(n), ev_keys_.begin());
    thrust::copy(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), ev_perm_.begin());
    thrust::stable_sort_by_key(w.on_this(), ev_keys_.begin(), ev_keys_.end(), ev_perm_.begin());
    const auto perm  = thrust::raw_pointer_cast(ev_perm_.data());   // event of sorted position i
    const auto keys_p = thrust::raw_pointer_cast(ev_keys_.data());   // its cell
    const auto qs    = thrust::make_transform_iterator(ord0, [=](const I i) { return q[perm[i]]; });
    // running sum in sorted order, then the counts written back to the events' own order
    thrust::inclusive_scan_by_key(
      w.on_this(), keys_p, keys_p + static_cast<std::ptrdiff_t>(n), qs, cs);
    thrust::for_each(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      const auto j    = perm[i];
      const auto base = car_p[keys_p[i]];
      const auto m    = sstd::floor(base + cs[i]) - sstd::floor(base + cs[i] - q[j]);
      out[j]          = static_cast<std::int32_t>(m);
    });
  }
  // one lane per event: the cell's remainder advances by what was asked minus what was given
  thrust::for_each(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
    sstd::atomic_add(&car_p[cell[i]], q[i] - static_cast<double>(out[i]));
  });
}

// ---------------------------------------------------------------- the step
// the padded layout of every species, each on its own random stream
void
  Sampler::layout(TileView& t)
{
  const auto lp = LayoutParams { .nbin       = par.nbin,
                                 .ncell_edge = par.ncell_edge,
                                 .slack      = par.slack,
                                 .smin       = par.smin };
  for(auto s = 0uz; s < 3uz; ++s) {
    t.pop[s]->layout(*t.con[s], lp, t.mins, t.maxs, runko::rng_key(par.seed, istep_, s, 0));
  }
}

// a layout when any species lost its tables (a source on a fresh or reset tile), or on the
// epoch beat of the caller
void
  Sampler::lay_if_needed(TileView& t, const bool epoch)
{
  auto need = epoch;
  for(auto s = 0uz; s < 3uz; ++s) { need = need or not t.pop[s]->laid(); }
  if(need) { this->layout(t); }
}

void
  Sampler::step(TileView& t)
{
  // ---- a new epoch every nsort steps, or as soon as a layout went stale
  const auto epoch = istep_ % static_cast<std::uint64_t>(par.nsort) == 0;
  this->lay_if_needed(t, epoch);
  if(not epoch) {
    // whatever moved the particles since the last step (the PIC push) did not widen the key
    // bounds; one sweep keeps every envelope a majorant (a no-op when nothing moved untouched)
    for(auto s = 0uz; s < 3uz; ++s) { t.pop[s]->touch_all(*t.con[s]); }
  }

  // ---- the context of the step and its slot/live statistics
  Ctx c { t, par.nbin };
  for(auto s = 0uz; s < 3uz; ++s) {
    stats.slots[s] += static_cast<double>(t.con[s]->size());
    stats.live[s] += static_cast<double>(t.con[s]->number_alive());
    t.ids_used[s] = 0;
  }
  thrust::fill(counters_.begin(), counters_.end(), 0);

  // The channels run one after another in a random order, each reading the states the
  // previous one left; the splitting error is second order in the per-channel rates and the
  // random order cancels it in expectation.
  auto order = std::vector<int> {};   // the enabled channels, shuffled below
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { order.push_back(static_cast<int>(i)); }
  }
  auto rng = runko::Rng { runko::rng_key(par.seed, istep_, 0, 0), 0 };
  // Fisher-Yates shuffle on the step's own stream
  for(auto i = order.size(); i > 1; --i) {
    const auto j = static_cast<std::size_t>(rng.integer(i));
    std::swap(order[i - 1], order[j]);
  }
  const auto b_tile = t.field.b_tile;   // the tile field bounds the absorption envelopes
  // dispatch each channel to the kernels of its record type
  for(const auto l: order) {
    switch(static_cast<Channel>(l)) {
      case Channel::compton_m: twobody(c, ChCompton { Species::em }, l); break;
      case Channel::compton_p: twobody(c, ChCompton { Species::ep }, l); break;
      case Channel::pairann: twobody(c, ChAnn {}, l); break;
      case Channel::photann: twobody(c, ChGG {}, l); break;
      case Channel::ssa_m: twobody(c, ChSSA { Species::em, b_tile, box.x_min }, l); break;
      case Channel::ssa_p: twobody(c, ChSSA { Species::ep, b_tile, box.x_min }, l); break;
      default: onebody(c, static_cast<Channel>(l), l); break;
    }
  }
  ++stats.steps;
  ++istep_;
}

// ---------------------------------------------------------------- two-body channel
template<typename Ch>
void
  Sampler::twobody(Ctx& c, const Ch& ch, const int l)
{
  // ---- the channel's species and the sizes of its tables
  const auto i1    = static_cast<std::size_t>(Ch::t1);   // incident species index
  const auto i2    = static_cast<std::size_t>(ch.t2);    // partner species index
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;                          // (cell, bin) groups
  const auto npair = ncb * nb;                            // (cell, bin a, bin b) pairs
  auto& st         = stats.ch[static_cast<std::size_t>(l)];
  const auto w     = tyvi::mdgrid_work {};
  c.refresh_all();

  // ---- the tables of the two sides: per (cell, bin) the sub-ranges and the bounds over them
  const auto side = [&](
                      const std::size_t is,
                      const Sel sel,
                      thrust::device_vector<I>& base0,
                      thrust::device_vector<I>& cnt0,
                      thrust::device_vector<float>& wmax,
                      thrust::device_vector<float>& emin,
                      thrust::device_vector<float>& emax) {
    base0.resize(static_cast<std::size_t>(2 * ncb));
    cnt0.resize(static_cast<std::size_t>(2 * ncb));
    wmax.resize(static_cast<std::size_t>(ncb));
    emin.resize(static_cast<std::size_t>(ncb));
    emax.resize(static_cast<std::size_t>(ncb));
    const auto lay = c.v[is].lay;
    const auto b0 = thrust::raw_pointer_cast(base0.data()), c0 = thrust::raw_pointer_cast(cnt0.data());  // [0, ncb) sub-range 0, [ncb, 2 ncb) sub-range 1
    const auto wm = thrust::raw_pointer_cast(wmax.data());
    const auto lo = thrust::raw_pointer_cast(emin.data()), hi = thrust::raw_pointer_cast(emax.data());
    const auto inf = std::numeric_limits<float>::infinity();
    // one lane per (cell, bin): the keys the side spans and the bounds over them
    thrust::for_each(
      w.on_this(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(ncb), [=](const int cb) {
        int k0 = -1, k1 = -1;   // the one or two keys of this group
        switch(sel) {
          case Sel::cell: k0 = cb; break;
          case Sel::bit0: k0 = 2 * cb; break;
          case Sel::bit1: k0 = 2 * cb + 1; break;
          case Sel::both: k0 = 2 * cb; k1 = 2 * cb + 1; break;
        }
        // fill records demand and exceeds cap when a key overflowed: only the fitted slots
        // are proposable (past cap lives the next key, or the tail)
        b0[cb] = lay.b[k0];
        c0[cb] = sstd::min(lay.fill[k0], lay.cap[k0]);
        auto wmx = lay.wmax[k0], emn = lay.emin[k0], emx = lay.emax[k0];
        if(k1 >= 0) {
          b0[ncb + cb] = lay.b[k1];
          c0[ncb + cb] = sstd::min(lay.fill[k1], lay.cap[k1]);
          wmx = sstd::max(wmx, lay.wmax[k1]);
          emn = sstd::min(emn, lay.emin[k1]);
          emx = sstd::max(emx, lay.emax[k1]);
        } else {
          b0[ncb + cb] = 0;
          c0[ncb + cb] = 0;
        }
        wm[cb] = wmx;
        lo[cb] = c0[cb] + c0[ncb + cb] > 0 ? emn : inf;
        hi[cb] = emx;
      });
  };
  side(i1, Ch::sel1, base1_, cnt1_, wmax1_, emin1_, emax1_);
  side(i2, Ch::sel2, base2_, cnt2_, wmax2_, emin2_, emax2_);
  // the two sides as the kernels read them
  const auto s1 = SideView { thrust::raw_pointer_cast(base1_.data()), thrust::raw_pointer_cast(cnt1_.data()),
                             thrust::raw_pointer_cast(base1_.data()) + ncb, thrust::raw_pointer_cast(cnt1_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax1_.data()), thrust::raw_pointer_cast(emin1_.data()),
                             thrust::raw_pointer_cast(emax1_.data()) };
  const auto s2 = SideView { thrust::raw_pointer_cast(base2_.data()), thrust::raw_pointer_cast(cnt2_.data()),
                             thrust::raw_pointer_cast(base2_.data()) + ncb, thrust::raw_pointer_cast(cnt2_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax2_.data()), thrust::raw_pointer_cast(emin2_.data()),
                             thrust::raw_pointer_cast(emax2_.data()) };

  // ---- the envelope of every bin pair: slot counts and the bounds of those bins, nothing else
  Lam_.resize(static_cast<std::size_t>(npair));
  sh_.resize(static_cast<std::size_t>(npair));
  wm_.resize(static_cast<std::size_t>(npair));
  rowsum_.resize(static_cast<std::size_t>(ncb));
  colsum_.resize(static_cast<std::size_t>(ncb));
  const auto Lam = thrust::raw_pointer_cast(Lam_.data()), sh = thrust::raw_pointer_cast(sh_.data());  // envelope and rate bound per bin pair
  const auto wmp = thrust::raw_pointer_cast(wm_.data());                                              // weight bound per bin pair
  const auto rowsum = thrust::raw_pointer_cast(rowsum_.data()), colsum = thrust::raw_pointer_cast(colsum_.data());  // Lam summed over b, over a
  const auto kappa = par.kappa;
  const auto N2c   = box.N2() / ncell;   // the two-body clock of one reaction cell
  // one lane per incident (cell, bin a): the envelope of every partner bin b and their sum
  thrust::for_each(
    w.on_this(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(ncb), [=](const int ca) {
      const auto cc  = ca / nb;
      const auto n1  = static_cast<double>(s1.N(ca));
      const auto lo1 = static_cast<double>(s1.emin[ca]), hi1 = static_cast<double>(s1.emax[ca]);
      const auto w1  = static_cast<double>(s1.wmax[ca]);
      auto sum       = 0.0;
      for(auto bb = 0; bb < nb; ++bb) {
        const auto cb  = cc * nb + bb;
        const auto g   = ca * nb + bb;
        const auto n2  = static_cast<double>(s2.N(cb));
        const auto lo2 = static_cast<double>(s2.emin[cb]), hi2 = static_cast<double>(s2.emax[cb]);
        const auto w2  = static_cast<double>(s2.wmax[cb]);
        auto shv       = n1 > 0.0 and n2 > 0.0 ? ch.bound(lo1, hi1, lo2, hi2) : 0.0;
        // bound on the effective weight of the merged pairs from the two maxima alone: f_kin is
        // extremal at the corners of the gamma range, so they bound it from below; leptons
        // below the channel's threshold never react, so the range starts there
        auto wmv = sstd::max(w1, w2);
        if constexpr(Ch::has_dgm) {
          if(kappa > 0.0) {
            const auto glo = sstd::max(lo2, Ch::g_min);
            const auto fk  = sstd::min(fkin(ch, kappa, hi1, glo), fkin(ch, kappa, hi1, hi2));
            wmv           = sstd::max(w2, w1 / fk);
          }
        }
        auto L = Ch::c_l * n1 * n2 * 2.0 * shv * wmv / N2c;
        if(not(L == L) or L > 1e300) { L = 0.0; }
        Lam[g] = L;
        sh[g]  = shv;
        wmp[g] = wmv;
        sum += L;
      }
      rowsum[ca] = sum;
    });
  // one lane per partner (cell, bin b): the envelope summed over the incident bins
  thrust::for_each(
    w.on_this(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(ncb), [=](const int cb) {
      const auto cc = cb / nb, bb = cb % nb;
      auto sum      = 0.0;
      for(auto a = 0; a < nb; ++a) { sum += Lam[(cc * nb + a) * nb + bb]; }
      colsum[cb] = sum;
    });
  stats.cell_lanes += npair;
  // the envelope-weighted mean of the rate bound, a measure of the bins' tightness
  const auto g0 = thrust::counting_iterator<int>(0);
  st.shbar_num += thrust::transform_reduce(
    w.on_this(), g0, g0 + npair, [=](const int g) { return Lam[g] * sh[g]; }, 0.0, thrust::plus<double> {});
  st.shbar_den += thrust::transform_reduce(
    w.on_this(), g0, g0 + npair, [=](const int g) { return Lam[g]; }, 0.0, thrust::plus<double> {});

  // Distinct strata need at most as many events per bin as the bin has slots, on both sides;
  // where the envelope exceeds that the channel is split into sub-steps of a larger N_2. The
  // sub-steps of the additive channels are capped: beyond the cap the row clamp drops the
  // surplus events (an under-count the cap controls), so no two lanes ever hold one slot.
  // sub-steps from the incident side: events of a bin over its slots less one
  auto nsub = thrust::transform_reduce(
    w.on_this(), g0, g0 + ncb,
    [=](const int ca) {
      const auto n = static_cast<double>(s1.N(ca));
      return rowsum[ca] > 0.0 and n >= 1.0 ? sstd::ceil(rowsum[ca] / sstd::max(n - 1.0, 1.0)) : 1.0;
    },
    1.0, thrust::maximum<double> {});
  if constexpr(Ch::kind == Kind::convert) {
    // a consumed partner needs distinct strata too: the same over the partner side
    nsub = sstd::max(nsub, thrust::transform_reduce(
      w.on_this(), g0, g0 + ncb,
      [=](const int cb) {
        const auto n = static_cast<double>(s2.N(cb));
        return colsum[cb] > 0.0 and n >= 1.0 ? sstd::ceil(colsum[cb] / sstd::max(n - 1.0, 1.0)) : 1.0;
      },
      1.0, thrust::maximum<double> {}));
  } else {
    nsub = sstd::min(nsub, static_cast<double>(par.nsub_max));
  }
  const auto ns = static_cast<int>(nsub);
  st.nsub_max_seen = sstd::max(st.nsub_max_seen, ns);
  st.nsub_steps += (ns > 1);
  // the sub-steps, each with the envelope divided by their number
  for(auto r = 0; r < ns; ++r) { pass(c, ch, l, r, ns, s1, s2); }
}

template<typename Ch>
void
  Sampler::pass(Ctx& c, const Ch& ch, const int l, const int r, const int nsub, const SideView& s1, const SideView& s2)
{
  // ---- sizes, one random stream per kernel of this (channel, sub-step), the device stats
  const auto i1    = static_cast<std::size_t>(Ch::t1);   // incident species index
  const auto i2    = static_cast<std::size_t>(ch.t2);    // partner species index
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;                          // (cell, bin) groups
  const auto npair = ncb * nb;                            // (cell, bin a, bin b) pairs
  auto& st         = stats.ch[static_cast<std::size_t>(l)];
  const auto w     = tyvi::mdgrid_work {};
  const auto step  = istep_;
  const auto keyA  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(1 + 4 * l));  // cell kernels: counts, shifts
  const auto keyB  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(2 + 4 * l));  // thinning kernel: strata, acceptance
  const auto keyC  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(3 + 4 * l));  // react kernels: products, copies
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ds = thrust::raw_pointer_cast(dstat_.data());   // the DStat accumulators
  const auto Lam = thrust::raw_pointer_cast(Lam_.data()), sh = thrust::raw_pointer_cast(sh_.data());  // envelope and rate bound per bin pair
  const auto wmp = thrust::raw_pointer_cast(wm_.data());                                              // weight bound per bin pair

  // ---- cell kernel: counts by one cumulative rounding per incident bin, clamps, shifts
  M_.resize(static_cast<std::size_t>(npair));
  off1_.resize(static_cast<std::size_t>(npair));
  off2_.resize(static_cast<std::size_t>(npair));
  kk_.resize(static_cast<std::size_t>(npair));
  M1_.resize(static_cast<std::size_t>(ncb));
  M2_.resize(static_cast<std::size_t>(ncb));
  sh1_.resize(static_cast<std::size_t>(ncb));
  sh2_.resize(static_cast<std::size_t>(ncb));
  const auto M = thrust::raw_pointer_cast(M_.data()), off1 = thrust::raw_pointer_cast(off1_.data());   // events per bin pair; stratum base along b
  const auto off2 = thrust::raw_pointer_cast(off2_.data()), kk = thrust::raw_pointer_cast(kk_.data());  // stratum base along a; partner stratum offset
  const auto M1 = thrust::raw_pointer_cast(M1_.data()), M2 = thrust::raw_pointer_cast(M2_.data());      // events per incident bin, per partner bin
  const auto sh1 = thrust::raw_pointer_cast(sh1_.data()), sh2 = thrust::raw_pointer_cast(sh2_.data());  // cyclic shift of each range
  const auto cb0 = thrust::counting_iterator<int>(0);
  // one lane per incident (cell, bin a): the counts M[a, b] by cumulative rounding, row clamp
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int ca) {
    auto rng  = runko::Rng { keyA, static_cast<std::uint64_t>(ca) };
    const auto U = rng.uniform();     // one uniform rounds the whole row
    auto C = 0.0, Fp = 0.0;           // cumulative envelope and its previous floor
    for(auto bb = 0; bb < nb; ++bb) {
      const auto g = ca * nb + bb;
      C += Lam[g] / nsub;
      const auto F = sstd::floor(C + U);
      M[g]         = static_cast<std::int32_t>(F - Fp);
      Fp           = F;
    }
    {
      // hold the events of the row below its slot count, so the strata stay distinct (a
      // backstop the sub-cycling makes rare); the drops are counted as clamped
      auto ex = 0;
      for(auto bb = 0; bb < nb; ++bb) { ex += M[ca * nb + bb]; }
      ex -= static_cast<int>(s1.N(ca));
      for(; ex > 0; --ex) {
        auto imax = 0;
        for(auto bb = 1; bb < nb; ++bb) {
          if(M[ca * nb + bb] > M[ca * nb + imax]) { imax = bb; }
        }
        --M[ca * nb + imax];
        sstd::atomic_add(&ds[d_clamped], 1.0);
      }
    }
  });
  // one lane per partner (cell, bin b): column clamp, strata bases along a, the partner shift
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int cb) {
    const auto cc = cb / nb, bb = cb % nb;
    if constexpr(Ch::kind == Kind::convert) {
      // hold the events of the column below the partner bin's slot count
      auto ex = 0;
      for(auto a = 0; a < nb; ++a) { ex += M[(cc * nb + a) * nb + bb]; }
      ex -= static_cast<int>(s2.N(cb));
      for(; ex > 0; --ex) {
        auto imax = 0;
        for(auto a = 1; a < nb; ++a) {
          if(M[(cc * nb + a) * nb + bb] > M[(cc * nb + imax) * nb + bb]) { imax = a; }
        }
        --M[(cc * nb + imax) * nb + bb];
        sstd::atomic_add(&ds[d_clamped], 1.0);
      }
    }
    auto sum = 0;
    for(auto a = 0; a < nb; ++a) {
      const auto g = (cc * nb + a) * nb + bb;
      off2[g]      = sum;  // strata bases: disjoint across bin pairs
      sum += M[g];
    }
    M2[cb]   = sum;
    auto rng = runko::Rng { keyA, static_cast<std::uint64_t>(ncb + cb) };
    sh2[cb]  = static_cast<std::int32_t>(rng.uniform() * static_cast<double>(sstd::max(s2.N(cb), I { 1 })));
  });
  // one lane per incident (cell, bin a): strata bases along b, stratum offsets, incident shift
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int ca) {
    auto rng = runko::Rng { keyA, static_cast<std::uint64_t>(2 * ncb + ca) };
    auto sum = 0;
    for(auto bb = 0; bb < nb; ++bb) {
      const auto g = ca * nb + bb;
      off1[g]      = sum;
      sum += M[g];
      kk[g] = static_cast<std::int32_t>(rng.uniform() * static_cast<double>(sstd::max(M[g], 1)));
    }
    M1[ca]  = sum;
    sh1[ca] = static_cast<std::int32_t>(rng.uniform() * static_cast<double>(sstd::max(s1.N(ca), I { 1 })));
  });
  // first lane of every bin pair from the running sum of the counts; the last entry is the
  // lane total of the pass, read on the host once the count kernels and the scan are done
  Mstart_.resize(static_cast<std::size_t>(npair) + 1);
  Mstart_[0] = I { 0 };
  thrust::inclusive_scan(w.on_this(), M_.begin(), M_.end(), Mstart_.begin() + 1);
  w.wait();
  const auto Mtot = static_cast<std::size_t>(Mstart_[static_cast<std::size_t>(npair)]);
  st.M += static_cast<double>(Mtot);
  if(Mtot == 0) { return; }

  // ---- thinning kernel: one lane per envelope event
  lane_acc_.resize(Mtot);
  lane_i_.resize(Mtot);
  lane_j_.resize(Mtot);
  lane_f_.resize(Mtot);
  const auto Mstart = thrust::raw_pointer_cast(Mstart_.data());
  const auto lacc = thrust::raw_pointer_cast(lane_acc_.data());                                      // accepted flag per lane
  const auto li = thrust::raw_pointer_cast(lane_i_.data()), lj = thrust::raw_pointer_cast(lane_j_.data());  // incident and partner slot per lane
  const auto lf = thrust::raw_pointer_cast(lane_f_.data());                                          // accumulation factor per lane
  const auto v1 = c.v[i1], v2 = c.v[i2];   // views of the incident and partner species
  const auto fv = c.t->field;              // the field at the partner, for the absorption rate
  const auto kappa = par.kappa;
  const auto lane0 = thrust::counting_iterator<I>(0);
  // one lane per envelope event: draw the pair, test the weight factor, then the rate
  thrust::for_each(w.on_this(), lane0, lane0 + static_cast<std::ptrdiff_t>(Mtot), [=](const I lane) {
    auto rng     = runko::Rng { keyB, lane };
    const auto g = segment_of(Mstart, npair, lane);          // the lane's bin pair
    const auto e = static_cast<std::uint64_t>(lane - Mstart[g]);   // its index within the pair
    const auto cc = g / (nb * nb), a = (g / nb) % nb, bb = g % nb;  // cell, incident bin, partner bin
    const auto cb1 = cc * nb + a, cb2 = cc * nb + bb;               // the two (cell, bin) groups
    const auto Mg = static_cast<std::uint64_t>(M[g]);               // events of the pair
    const auto n1 = s1.N(cb1), n2 = s2.N(cb2);                       // slots of the two groups
    const auto p1 = static_cast<std::uint64_t>(off1[g]) + e;        // incident stratum
    const auto vv1 = (static_cast<I>(sh1[cb1]) + stratum(p1, static_cast<std::uint64_t>(M1[cb1]), n1, rng)) % n1;  // virtual index on side 1
    I vv2;                                                          // virtual index on side 2
    if constexpr(Ch::kind == Kind::convert) {
      // a consumed partner: distinct strata, offset within the pair to decorrelate from p1
      const auto p2 = static_cast<std::uint64_t>(off2[g]) + (e + static_cast<std::uint64_t>(kk[g])) % Mg;
      vv2 = (static_cast<I>(sh2[cb2]) + stratum(p2, static_cast<std::uint64_t>(M2[cb2]), n2, rng)) % n2;
    } else {
      // the partner only accumulates what the events give it, so it may be the partner of
      // several lanes: drawn uniformly with replacement
      vv2 = static_cast<I>(rng.uniform() * static_cast<double>(n2)) % n2;
    }
    const auto i = s1.pick(cb1, vv1), j = s2.pick(cb2, vv2);   // the two slots
    li[lane] = i;
    lj[lane] = j;
    lf[lane] = 1.0;
    lacc[lane] = 0;
    // a dead slot, or one dropped from its key by a kick: nothing to react
    const auto virt = not(v1.alive(i) and v2.alive(j)) or v1.lay.skey[i] < 0 or v2.lay.skey[j] < 0;
    if(virt) {
      sstd::atomic_add(&ds[d_virt], 1.0);
      return;
    }
    // Two acceptance stages against one uniform: the weight factor first, at the largest value
    // the merged pair can have; then the effective weight of the merged pair and the cross
    // section, whose threshold is the first times a factor <= 1.
    const auto u  = rng.uniform();                        // the one uniform of both stages
    const auto wi = v1.w(i), wj = v2.w(j);
    const auto wm = wmp[g];                               // weight bound of the bin pair
    const auto wf = sstd::min(sstd::max(wi, wj) / wm, 1.0);   // weight-stage threshold
    if(not(u < wf)) { return; }
    sstd::atomic_add(&ds[d_wpass], 1.0);
    const auto p1v = v1.u(i), p2v = v2.u(j);
    const auto rate = pair_rate(ch, fv, p1v, p2v);
    auto weff = sstd::max(wi, wj);                        // effective weight of the merged pair
    if constexpr(Ch::has_dgm) {
      if(kappa > 0.0) {
        // heavy-light merging: f events at once, the effective weight max(w_j, w_i/f_kin)
        const auto x   = toolbox::norm(p1v);
        const auto gam = sstd::sqrt(1.0 + toolbox::dot(p2v, p2v));
        weff           = sstd::max(wj, wi / fkin(ch, kappa, x, gam));
        lf[lane]       = sstd::max(wi, wj) / weff;
      }
    }
    if(u < (weff / wm) * rate / sh[g]) { lacc[lane] = 1; }   // rate stage against the bound
  });
  // compact the accepted lanes into acc_
  acc_.resize(Mtot);
  const auto acc_end = thrust::copy_if(
    w.on_this(), lane0, lane0 + static_cast<std::ptrdiff_t>(Mtot), acc_.begin(), [=](const I lane) { return lacc[lane] != 0; });
  const auto na = static_cast<std::size_t>(acc_end - acc_.begin());   // accepted events
  st.acc += static_cast<double>(na);
  {
    // read back the lane statistics of the pass
    thrust::host_vector<double> h = dstat_;
    st.virt += h[d_virt];
    st.wpass += h[d_wpass];
    st.clamped += h[d_clamped];
  }
  if(na == 0) {
    w.wait();
    return;
  }

  // ---- react kernel, first half: the products of every accepted event and its decisions
  ev_q1_.resize(na);
  ev_q2_.resize(na);
  ev_q3_.resize(na);
  ev_n1_.resize(na);
  ev_n2_.resize(na);
  ev_n3_.resize(na);
  ev_cell_.resize(na);
  ev_p_.resize(6 * na);
  const auto acc = thrust::raw_pointer_cast(acc_.data());                                            // lane of every accepted event
  const auto q1 = thrust::raw_pointer_cast(ev_q1_.data()), q2 = thrust::raw_pointer_cast(ev_q2_.data());  // update/removal probability of parent 1, of parent 2
  const auto q3 = thrust::raw_pointer_cast(ev_q3_.data());                                           // lepton pairs to create (photon-photon)
  const auto n1p = thrust::raw_pointer_cast(ev_n1_.data()), n2p = thrust::raw_pointer_cast(ev_n2_.data());  // their carry-rounded decisions
  const auto n3p = thrust::raw_pointer_cast(ev_n3_.data());
  const auto ecell = thrust::raw_pointer_cast(ev_cell_.data());                                      // reaction cell of every event
  const auto ep = thrust::raw_pointer_cast(ev_p_.data());                                            // the two product momenta, 6 floats per event
  const auto wl = box.w_lep;
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ev0 = thrust::counting_iterator<I>(0);
  // one lane per accepted event: sample the products, set the decision probabilities
  thrust::for_each(w.on_this(), ev0, ev0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
    auto rng       = runko::Rng { keyC, e };
    const auto lane = acc[e];
    const auto i = li[lane], j = lj[lane];   // the two parents
    const auto f = lf[lane];                 // events merged into this one
    const auto g = segment_of(Mstart, npair, lane);
    ecell[e]     = g / (nb * nb);
    const auto p1v = v1.u(i), p2v = v2.u(j);
    const auto wi = v1.w(i), wj = v2.w(j);
    const auto s  = sstd::min(wi, wj);  // the weight one event transfers
    sstd::atomic_add(&ds[d_wmoved], f * s);
    sstd::atomic_add(&ds[d_facc], f);
    auto ntry = 0;        // proposals of the rejection sampler
    Vec3 pa, pb;          // the two products: (k1, kick), (f x, -), (k0, k1) or (e-, e+)
    if constexpr(Ch::kind == Kind::scatter) {
      // Compton: the photon takes the scattered state with the weight f s; the lepton takes f
      // exchanges in one, a momentum kick added to whatever other lanes give it
      const auto P = comp::sample(p1v, p2v, rng, ntry);
      pa           = P.k1;
      pb           = (-f) * (P.k1 - p1v);
      q1[e]        = sstd::min(f * s / wi, 1.0);  // the photon moves the merged weight f s
      q2[e]        = sstd::min(s / wj, 1.0);
      q3[e]        = 0.0;
    } else if constexpr(Ch::kind == Kind::absorb) {
      // self-absorption: the photon loses the weight f s, the lepton gains the energy f x
      pa    = Vec3(f * toolbox::norm(p1v), 0.0, 0.0);
      pb    = Vec3(0.0, 0.0, 0.0);
      q1[e] = sstd::min(f * s / wi, 1.0);  // the photon loses the merged weight f s
      q2[e] = sstd::min(s / wj, 1.0);
      q3[e] = 0.0;
    } else if constexpr(Ch::t2 == Species::ep) {
      // pair annihilation: two photons of weight s, both parents removed
      const auto P = ann::sample(p1v, p2v, rng, ntry);
      pa           = P.k0;
      pb           = P.k1;
      q1[e]        = sstd::min(s / wi, 1.0);
      q2[e]        = sstd::min(s / wj, 1.0);
      q3[e]        = 0.0;
    } else {
      // photon-photon pair production: the weight s becomes leptons at the lepton weight
      const auto P = gg::sample(p1v, p2v, rng, ntry);
      const auto flip = rng.uniform() < 0.5;
      pa              = flip ? P.p1 : P.p2;  // e-
      pb              = flip ? P.p2 : P.p1;  // e+
      q1[e]           = sstd::min(s / wi, 1.0);
      q2[e]           = sstd::min(s / wj, 1.0);
      q3[e]           = s / wl;  // pairs at the lepton weight, rounded systematically
    }
    // store the products for the second half
    for(auto k = 0uz; k < 3uz; ++k) {
      ep[6 * e + k]     = static_cast<float>(pa[k]);
      ep[6 * e + 3 + k] = static_cast<float>(pb[k]);
    }
    sstd::atomic_add(&ds[d_trials], static_cast<double>(ntry));
    sstd::atomic_add(&ds[d_trials_n], 1.0);
  });
  w.wait();
  {
    // read back the event statistics
    thrust::host_vector<double> h = dstat_;
    st.wmoved += h[d_wmoved];
    st.facc += h[d_facc];
    st.trials += h[d_trials];
    st.trials_n += h[d_trials_n];
  }

  // ---- the decisions, by the systematic rounding of every cell
  const auto tag1 = std::to_string(static_cast<int>(Ch::t1)), tag2 = std::to_string(static_cast<int>(ch.t2));   // carry sites are per species
  // scatter: update the photon, kick the lepton; absorb: remove the photon, kick the lepton;
  // convert: remove both, and for photon-photon also count the lepton pairs
  if constexpr(Ch::kind == Kind::scatter) {
    select(w, true, "sc" + tag1, ncell, na, ecell, q1, n1p);
    select(w, true, "kk" + tag2, ncell, na, ecell, q2, n2p);
  } else if constexpr(Ch::kind == Kind::absorb) {
    select(w, true, "rm" + tag1, ncell, na, ecell, q1, n1p);
    select(w, true, "kk" + tag2, ncell, na, ecell, q2, n2p);
  } else {
    select(w, true, "rm" + tag1, ncell, na, ecell, q1, n1p);
    select(w, true, "rm" + tag2, ncell, na, ecell, q2, n2p);
    if constexpr(Ch::t2 == Species::ph) { select(w, true, "lep", ncell, na, ecell, q3, n3p); }
  }

  // ---- react kernel, second half: parents updated in place, products staged
  staged_.resize(2 * na);                                     // at most two products per event
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());   // the Counter accumulators
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  touched_.assign(c.t->con[i2]->size(), std::uint8_t { 0 });   // partner slots kicked this pass
  const auto touched = thrust::raw_pointer_cast(touched_.data());
  if constexpr(Ch::kind == Kind::absorb) { dgam_.assign(c.t->con[i2]->size(), 0.0f); }   // absorbed energy per partner slot
  const auto dgam = thrust::raw_pointer_cast(dgam_.data());
  const auto a_res = par.a_res, x_ref = par.x_ref;
  const auto n_cap = par.n_cap;
  // append one product record to the staging buffer
  const auto stage = [=](const Vec3& p, const float* pos, const double w_each, const int n_r, const int cell, const int species, const int bit) {
    const auto k = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    stg[k]       = staged_at(p, pos, w_each, n_r, cell, species, bit);
  };
  // one lane per accepted event: apply the decisions to the parents, stage the products
  thrust::for_each(w.on_this(), ev0, ev0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
    auto rng        = runko::Rng { keyC, e + (std::uint64_t { 1 } << 32) };   // past the first half's lanes
    const auto lane = acc[e];
    const auto i = li[lane], j = lj[lane];
    const auto cell = ecell[e];
    const auto pa   = Vec3(static_cast<double>(ep[6 * e]), static_cast<double>(ep[6 * e + 1]), static_cast<double>(ep[6 * e + 2]));      // product 1
    const auto pb   = Vec3(static_cast<double>(ep[6 * e + 3]), static_cast<double>(ep[6 * e + 4]), static_cast<double>(ep[6 * e + 5]));  // product 2
    const float posi[3] = { v1.pos[i][0], v1.pos[i][1], v1.pos[i][2] };   // products are born where the incident is
    if constexpr(Ch::kind == Kind::scatter) {
      if(n1p[e] > 0) {
        // the parent takes the new state and is realized as the copies the resolution function
        // asks for; the parent is the first copy
        const auto x0 = v1.e(i), x1 = toolbox::norm(pa);
        const auto n  = f_res(x1, a_res, x_ref) / f_res(x0, a_res, x_ref);
        const auto cp = copies(n, v1.w(i), rng, n_cap);
        if(cp.n_r >= 1) {
          v1.set_u(i, pa);
          v1.set_w(i, cp.w_each);
          v1.touch(i);
          if(cp.n_r > 1) { stage(pa, posi, cp.w_each, cp.n_r - 1, cell, 0, v1.lay.skey[i] >= 0 ? v1.lay.skey[i] % 2 : 0); }
        } else {
          v1.kill(i);  // the roulette merged it away
        }
      }
      if(n2p[e] > 0) {
        // the lepton's kick, summed with those of the other lanes that hit it
        for(auto k = 0uz; k < 3uz; ++k) { sstd::atomic_add(&v2.vel[j][k], static_cast<float>(pb[k])); }
        touched[j] = 1;
      }
    } else if constexpr(Ch::kind == Kind::absorb) {
      if(n2p[e] > 0) {
        // the absorbed energy, summed per lepton and applied after the kernel
        sstd::atomic_add(&dgam[j], static_cast<float>(pa[0]));
        touched[j] = 1;
      }
      if(n1p[e] > 0) { v1.kill(i); }
    } else if constexpr(Ch::t2 == Species::ep) {
      // annihilation: each photon realized as the copies its resolution asks for
      const auto s  = sstd::min(v1.w(i), v2.w(j));
      const auto x0 = toolbox::norm(pa), x1 = toolbox::norm(pb);
      const auto ca = copies(f_res(x0, a_res, x_ref), s, rng, n_cap);
      if(ca.n_r >= 1) { stage(pa, posi, ca.w_each, ca.n_r, cell, 0, static_cast<int>(rng.uniform() < 0.5)); }
      const auto cb = copies(f_res(x1, a_res, x_ref), s, rng, n_cap);
      if(cb.n_r >= 1) { stage(pb, posi, cb.w_each, cb.n_r, cell, 0, static_cast<int>(rng.uniform() < 0.5)); }
      if(n1p[e] > 0) { v1.kill(i); }
      if(n2p[e] > 0) { v2.kill(j); }
    } else {
      const auto nl = n3p[e];
      if(nl > 0) {
        // both leptons at the position of the incident photon: a pair born on top of itself adds
        // no charge density, so div E stays what the fields had
        stage(pa, posi, wl, nl, cell, 1, 0);
        stage(pb, posi, wl, nl, cell, 2, 0);
      }
      if(n1p[e] > 0) { v1.kill(i); }
      if(n2p[e] > 0) { v2.kill(j); }
    }
  });
  w.wait();
  // ---- the kicked leptons take their summed changes, then the products are placed
  if constexpr(Ch::kind == Kind::absorb) { apply_dgam(c, ch.t2); }
  if constexpr(Ch::kind == Kind::scatter) { touch_marked(c, ch.t2); }
  this->place(c);
}

// the lepton takes the summed energy of its absorptions along its direction of motion
void
  Sampler::apply_dgam(Ctx& c, const Species s)
{
  const auto is = static_cast<std::size_t>(s);
  const auto v  = c.v[is];
  const auto dgam = thrust::raw_pointer_cast(dgam_.data());
  const auto touched = thrust::raw_pointer_cast(touched_.data());
  // one lane per slot: rescale |u| to the energy after the absorptions, on the mass shell
  tyvi::mdgrid_work {}
    .for_each_index(v.ids, [=](const auto idx) {
      const auto i = static_cast<I>(idx[0]);
      if(touched[i] == 0 or not v.alive(i)) { return; }
      const auto dg = static_cast<double>(dgam[i]);
      if(dg <= 0.0) { return; }
      const auto u0 = v.u(i);
      const auto g0 = sstd::sqrt(1.0 + toolbox::dot(u0, u0));
      const auto g1 = g0 + dg;
      v.set_u(i, u0 * (sstd::sqrt(g1 * g1 - 1.0) / sstd::sqrt(sstd::max(g0 * g0 - 1.0, 1e-300))));
      v.touch(i);
    })
    .wait();
}

// widen the key bounds of the slots the kicks moved, once all kicks are in
void
  Sampler::touch_marked(Ctx& c, const Species s)
{
  const auto is = static_cast<std::size_t>(s);
  const auto v  = c.v[is];
  const auto touched = thrust::raw_pointer_cast(touched_.data());
  // one lane per slot
  tyvi::mdgrid_work {}
    .for_each_index(v.ids, [=](const auto idx) {
      const auto i = static_cast<I>(idx[0]);
      if(touched[i] != 0 and v.alive(i)) { v.touch(i); }
    })
    .wait();
}

// ---------------------------------------------------------------- products: slack of their key, overflow to the tail
void
  Sampler::place(Ctx& c)
{
  const auto nst = static_cast<std::size_t>(counters_[c_staged]);   // staged records
  if(nst == 0) { return; }
  const auto w   = tyvi::mdgrid_work {};
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto r0  = thrust::counting_iterator<I>(0);
  // copies per species, then room for them at the end of every container
  std::array<std::size_t, 3> tot {};   // copies to create per species
  for(auto s = 0; s < 3; ++s) {
    tot[static_cast<std::size_t>(s)] = static_cast<std::size_t>(thrust::transform_reduce(
      w.on_this(), r0, r0 + static_cast<std::ptrdiff_t>(nst),
      [=](const I k) { return stg[k].species == s ? static_cast<long>(stg[k].n_r) : 0L; }, 0L, thrust::plus<long> {}));
  }
  std::array<I, 3> tail_base {};   // first tail slot of every species: the old container end
  for(auto s = 0uz; s < 3uz; ++s) {
    tail_base[s] = static_cast<I>(c.t->con[s]->size());
    if(tot[s] > 0) {
      // grow once by the exact copy total; the layout follows with tail keys
      c.t->con[s]->grow(tot[s]);
      c.t->pop[s]->follow_container(*c.t->con[s]);
    }
  }
  c.refresh_all();   // the containers were resized
  thrust::fill(counters_.begin() + c_tail0, counters_.end(), 0);
  const auto cnt = thrust::raw_pointer_cast(counters_.data());   // tail and id counters
  const auto v   = c.v;
  const auto idb = c.t->id_base;                                  // first free id per species
  // one lane per staged record: every copy takes a slot in the slack of its key or in the tail
  thrust::for_each(w.on_this(), r0, r0 + static_cast<std::ptrdiff_t>(nst), [=](const I k) {
    const auto rec = stg[k];
    const auto s   = static_cast<std::size_t>(rec.species);
    const auto pv  = v[s];
    const auto lay = pv.lay;
    const auto e   = sstd::sqrt(static_cast<double>(pv.mass) * pv.mass + static_cast<double>(rec.u[0]) * rec.u[0] +
                              static_cast<double>(rec.u[1]) * rec.u[1] + static_cast<double>(rec.u[2]) * rec.u[2]);
    const auto kk  = lay.key_of(rec.cell, lay.bin_of(e), rec.bit);  // the product's own key
    for(auto n = 0; n < rec.n_r; ++n) {
      const auto f = sstd::atomic_fetch_add(&lay.fill[kk], I { 1 });   // rank behind the key's fill
      I slot;
      if(f < lay.cap[kk]) {
        slot = lay.b[kk] + f;  // behind what the key already holds, reactive at once
        // a product outside the bounds of its key must widen them, or the envelopes of the
        // channels still to come this step would no longer be upper bounds
        sstd::atomic_max(&lay.wmax[kk], rec.w);
        sstd::atomic_max(&lay.emax[kk], static_cast<float>(e));
        sstd::atomic_min(&lay.emin[kk], static_cast<float>(e));
        sstd::atomic_add(&cnt[c_slack], std::int32_t { 1 });
      } else {
        const auto t = sstd::atomic_fetch_add(&cnt[c_tail0 + rec.species], std::int32_t { 1 });
        slot         = tail_base[s] + static_cast<I>(t);  // inert until the next sort
        sstd::atomic_add(&cnt[c_tailn], std::int32_t { 1 });
      }
      sstd::atomic_add(&lay.prod[kk], I { 1 });  // demand, whether or not it fitted
      const auto id = idb[s] + static_cast<runko::prtc_id_type>(sstd::atomic_fetch_add(&cnt[c_id0 + rec.species], std::int32_t { 1 }));
      // write the copy
      for(auto q = 0uz; q < 3uz; ++q) {
        pv.pos[slot][q] = rec.pos[q];
        pv.vel[slot][q] = rec.u[q];
      }
      pv.ids[slot][] = id;
      if(pv.weighted) { pv.wgt[slot][] = rec.w; }
    }
  });
  w.wait();
  {
    // ids consumed and placement statistics
    thrust::host_vector<std::int32_t> h = counters_;
    for(auto s = 0uz; s < 3uz; ++s) {
      c.t->ids_used[s] += static_cast<std::size_t>(h[c_id0 + static_cast<int>(s)]);
      c.t->id_base[s] += static_cast<runko::prtc_id_type>(h[c_id0 + static_cast<int>(s)]);
    }
    stats.slackprod += h[c_slack];
    stats.tailprod += h[c_tailn];
  }
}

// ---------------------------------------------------------------- one-body channels
void
  Sampler::onebody(Ctx& c, const Channel which, const int l)
{
  // ---- the channel's species, the field and the one-body clock
  const auto emit = which != Channel::bw;                 // emission (syn) or pair production (bw)
  const auto is   = static_cast<std::size_t>(which == Channel::syn_m ? Species::em : which == Channel::syn_p ? Species::ep : Species::ph);
  auto& st        = stats.ch[static_cast<std::size_t>(l)];
  const auto w    = tyvi::mdgrid_work {};
  c.refresh_all();
  const auto v    = c.v[is];
  const auto n    = c.t->con[is]->size();                 // one lane per slot
  if(n == 0) { return; }
  const auto fv   = c.t->field;                           // the field at every particle
  const auto N1   = box.N1();
  const auto x_min = box.x_min;
  if(not emit) {
    // no photon above the pair threshold: nothing to fire
    const auto lay = v.lay;
    const auto k0  = thrust::counting_iterator<int>(0);
    const auto emax = thrust::transform_reduce(
      w.on_this(), k0, k0 + lay.nkeys, [=](const int k) { return lay.emax[k]; }, 0.0f, thrust::maximum<float> {});
    if(emax < static_cast<float>(bw::XMIN)) { return; }
  }
  // ---- random streams, device stats and the staging buffer of the pass
  const auto keyA = runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(100 + 2 * l));   // firing
  const auto keyC = runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(101 + 2 * l));   // products, copies
  st.M += static_cast<double>(n);
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ds = thrust::raw_pointer_cast(dstat_.data());
  staged_.resize(2 * n);                                    // at most two products per slot
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  const auto a_res = par.a_res, x_ref = par.x_ref, acc_syn = par.acc_syn, kappa = par.kappa;
  const auto n_cap = par.n_cap;
  const auto wl    = box.w_lep;
  // append one product record at the parent's position and cell
  const auto stage = [=](const Vec3& p, const I i, const double w_each, const int n_r, const int species, const int bit) {
    const auto k = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    const float pos[3] = { v.pos[i][0], v.pos[i][1], v.pos[i][2] };
    stg[k]       = staged_at(p, pos, w_each, n_r, v.lay.cell_of(pos[0], pos[1], pos[2]), species, bit);
  };
  const auto s0 = thrust::counting_iterator<I>(0);
  if(emit) {
    // Below one emission per step the event is resolved: it fires with the merged probability and
    // emits one photon of the merged weight. Above it the parent cools continuously every step and
    // emits one photon LP per step standing for all the emissions of the step; the resolution
    // function then thins the soft emission by the roulette of copies().
    // one lane per lepton slot: fire, cool, emit one photon LP
    thrust::for_each(w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      if(not v.alive(i)) {
        sstd::atomic_add(&ds[d_virt], 1.0);
        return;
      }
      auto rngA     = runko::Rng { keyA, i };
      auto rngC     = runko::Rng { keyC, i };
      const auto u  = v.u(i);
      const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
      const auto b   = fv.b_eff(v.pos[i][0], v.pos[i][1], v.pos[i][2], u / gam);   // the lepton's field
      const auto lam = syn::rate(u, b, x_min, false) / N1;   // emissions per step above the floor
      if(lam <= 0.0) { return; }
      const auto xm  = syn::xmean(gam * b, gam, false, x_min);   // mean photon energy above the floor
      const auto z   = sstd::max(rngA.uniform(), 1e-300);
      const auto slow = lam <= 1.0;                              // resolved: below one emission per step
      const auto f    = slow ? facc_syn(gam, xm, acc_syn, kappa) : 1.0;   // emissions merged into one LP
      const auto fire = slow ? -sstd::log(z) < lam / f : true;
      auto dgu = 0.0;                                            // energy the continuous cooling took
      if(not slow) {
        // the unresolved parent cools whether or not it emits an LP this step
        const auto g1 = sstd::max(gamma_cooled(gam, lam, xm, b, N1), 1.0);
        v.set_u(i, u * (sstd::sqrt(sstd::max(g1 * g1 - 1.0, 0.0)) / sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
        dgu = gam - g1;  // the energy the step took, to be carried by the LP
      }
      if(not fire) {
        if(not slow) { v.touch(i); }
        return;
      }
      const auto mult = slow ? 1.0 : dgu / sstd::max(xm, 1e-300);   // emissions the LP stands for
      const auto P    = syn::sample(u, b, rngC, false, x_min);      // the photon, from the pre-cooling state
      const auto x    = toolbox::norm(P.k);
      // a resolved parent gives up exactly what its photon carries, clipped at the mass shell; an
      // unresolved one has already cooled, and its LP carries the merged loss of f steps
      const auto dg  = slow ? f * x : dgu;                           // energy the LP should carry
      const auto dgc = sstd::min(dg, gam - 1.0);                     // clipped at the mass shell
      const auto wem = v.w(i) * f * mult * (dg > 0.0 ? dgc / sstd::max(dg, 1e-300) : 1.0);   // emitted weight
      if(slow) {
        // a resolved parent recoils by what its photon carries
        const auto g1 = sstd::max(gam - dgc, 1.0);
        v.set_u(i, u * (sstd::sqrt(sstd::max(g1 * g1 - 1.0, 0.0)) / sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
      }
      v.touch(i);
      sstd::atomic_add(&ds[d_wpass], 1.0);  // fired
      sstd::atomic_add(&ds[d_wem], wem);
      sstd::atomic_add(&ds[d_facc], f);
      // the photon LP, split or rouletted by the resolution function
      const auto cp = copies(f_res(x, a_res, x_ref), wem, rngC, n_cap);
      if(cp.n_r >= 1) { stage(P.k, i, cp.w_each, cp.n_r, 0, static_cast<int>(rngC.uniform() < 0.5)); }
    });
    w.wait();
    // read back the emission statistics
    thrust::host_vector<double> h = dstat_;
    st.virt += h[d_virt];
    st.acc += h[d_wpass];
    st.wem += h[d_wem];
    st.facc += h[d_facc];
  } else {
    // magnetic pair production: fire, then the pair at the lepton weight, rounded systematically
    ev_q3_.resize(n);          // pairs to create per fired photon, w/w_lep
    ev_n3_.resize(n);          // its carry-rounded count
    ev_cell_.resize(n);        // reaction cell of the fired photon
    lane_acc_.assign(n, 0);    // fired flag per slot
    const auto lacc = thrust::raw_pointer_cast(lane_acc_.data());
    // one lane per photon slot: fire on the conversion rate
    thrust::for_each(w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      if(not v.alive(i)) {
        sstd::atomic_add(&ds[d_virt], 1.0);
        return;
      }
      auto rngA      = runko::Rng { keyA, i };
      const auto k   = v.u(i);
      const auto b   = fv.b_eff(v.pos[i][0], v.pos[i][1], v.pos[i][2], k / sstd::max(toolbox::norm(k), 1e-300));   // the photon's field
      const auto lam = bw::rate(k, b) / N1;   // conversions per step
      const auto z   = sstd::max(rngA.uniform(), 1e-300);
      if(lam > 0.0 and -sstd::log(z) < lam) { lacc[i] = 1; }
    });
    // compact the fired slots into acc_
    acc_.resize(n);
    const auto acc_end = thrust::copy_if(
      w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), acc_.begin(), [=](const I i) { return lacc[i] != 0; });
    const auto na = static_cast<std::size_t>(acc_end - acc_.begin());   // fired photons
    {
      thrust::host_vector<double> h = dstat_;
      st.virt += h[d_virt];
    }
    st.acc += static_cast<double>(na);
    if(na == 0) {
      w.wait();
      return;
    }
    const auto acc = thrust::raw_pointer_cast(acc_.data());
    const auto q3  = thrust::raw_pointer_cast(ev_q3_.data());
    const auto n3  = thrust::raw_pointer_cast(ev_n3_.data());
    const auto ecell = thrust::raw_pointer_cast(ev_cell_.data());
    const auto e0 = thrust::counting_iterator<I>(0);
    // one lane per fired photon: the pair count it asks for and its cell
    thrust::for_each(w.on_this(), e0, e0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
      const auto i = acc[e];
      q3[e]        = v.w(i) / wl;
      ecell[e]     = v.lay.cell_of(v.pos[i][0], v.pos[i][1], v.pos[i][2]);
    });
    w.wait();
    // the counts by the systematic rounding of every cell
    select(w, false, "lep", c.ncell, na, ecell, q3, n3);
    // one lane per fired photon: draw the pair, stage its copies, remove the photon
    thrust::for_each(w.on_this(), e0, e0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
      const auto i = acc[e];
      auto rngC    = runko::Rng { keyC, i };
      const auto k = v.u(i);
      const auto b = fv.b_eff(v.pos[i][0], v.pos[i][1], v.pos[i][2], k / sstd::max(toolbox::norm(k), 1e-300));
      const auto P = bw::sample(k, b, rngC);
      if(n3[e] > 0) {
        // the pair at the photon's position, on top of each other: no charge density is added
        stage(P.um, i, wl, n3[e], 1, 0);
        stage(P.up, i, wl, n3[e], 2, 0);
      }
      v.kill(i);  // s = w: the parent always goes
    });
    w.wait();
  }
  this->place(c);
}

// ---------------------------------------------------------------- sources and sinks
namespace {

// slots of the per-call budget accumulated on the device (dbudget_): energy and weight
enum BStat : int { b_E = 0, b_W, b_count };

}  // namespace

void
  Sampler::inject_pairs_plaw(TileView& t, const double ell, const double p, const double gmin, const double gmax)
{
  this->lay_if_needed(t);
  // The number of pairs follows from the energy the compactness injects and the mean energy of the
  // spectrum, rounded stochastically; the two charges draw their energies independently.
  double gm;   // mean gamma of the power law
  if(p == 1.0) {
    gm = (gmax - gmin) / std::log(gmax / gmin);
  } else {
    const auto a = 1.0 - p, b = 2.0 - p;
    gm = p == 2.0 ? (b / a) * (std::pow(gmax, a) - std::pow(gmin, a)) / (std::pow(gmax, b) - std::pow(gmin, b))
                  : (a / b) * (std::pow(gmax, b) - std::pow(gmin, b)) / (std::pow(gmax, a) - std::pow(gmin, a));
  }
  const auto wl    = box.w_lep;
  const auto npair = box.dE(ell) / (2.0 * gm * wl);        // pairs per step at the lepton weight
  auto rng0        = runko::Rng { runko::rng_key(par.seed, istep_, 0, 300), 0 };
  const auto n     = static_cast<std::size_t>(std::floor(npair + rng0.uniform()));   // rounded stochastically
  if(n == 0) { return; }
  // ---- the staging buffer: one record per lepton, placed by the common place kernel
  Ctx c { t, par.nbin };
  staged_.resize(2 * n);
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  thrust::fill(dbudget_.begin(), dbudget_.end(), 0.0);
  const auto stg  = thrust::raw_pointer_cast(staged_.data());
  const auto cnt  = thrust::raw_pointer_cast(counters_.data());
  const auto db   = thrust::raw_pointer_cast(dbudget_.data());   // injected energy
  const auto key  = runko::rng_key(par.seed, istep_, 0, 301);
  const auto mins = t.mins, maxs = t.maxs;                       // tile extent for the positions
  const auto lay  = c.v[1].lay;
  const auto w    = tyvi::mdgrid_work {};
  const auto i0   = thrust::counting_iterator<I>(0);
  const auto glo  = sstd::pow(gmin, 1.0 - p), ghi = sstd::pow(gmax, 1.0 - p);   // the cumulative of gamma^-p at the ends
  // one lane per pair: a uniform position, the two energies from the power law, isotropic
  thrust::for_each(w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
    auto rng = runko::Rng { key, i };
    const float pos[3] = { mins[0] + (maxs[0] - mins[0]) * static_cast<float>(rng.uniform()),
                           mins[1] + (maxs[1] - mins[1]) * static_cast<float>(rng.uniform()),
                           mins[2] + (maxs[2] - mins[2]) * static_cast<float>(rng.uniform()) };
    const auto cell    = lay.cell_of(pos[0], pos[1], pos[2]);
    for(auto sp = 1; sp <= 2; ++sp) {  // electron and positron at the same position: neutral
      const auto z   = rng.uniform();
      const auto gam = sstd::pow(glo + z * (ghi - glo), 1.0 / (1.0 - p));   // inverse cumulative of gamma^-p
      const auto u   = sstd::sqrt(gam * gam - 1.0) * iso_dir(rng);
      const auto k   = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
      stg[k]         = staged_at(u, pos, wl, 1, cell, sp, 0);
      sstd::atomic_add(&db[b_E], wl * gam);
    }
  });
  w.wait();
  {
    // the energy and pair depth injected
    thrust::host_vector<double> h = dbudget_;
    budget.inj_E += h[b_E];
    budget.inj_W += 2.0 * static_cast<double>(n) * wl;
  }
  this->place(c);
}

void
  Sampler::inject_photons(TileView& t, const double ell, const double theta, const double wph, const bool planck)
{
  this->lay_if_needed(t);
  const auto xmean = (planck ? 2.701 : 3.0) * theta;       // mean photon energy of the spectrum
  const auto nph   = box.dE(ell) / (xmean * wph);          // photons per step at the weight wph
  auto rng0        = runko::Rng { runko::rng_key(par.seed, istep_, 0, 302), 0 };
  const auto n     = static_cast<std::size_t>(std::floor(nph + rng0.uniform()));   // rounded stochastically
  if(n == 0) { return; }
  // ---- the staging buffer: one record per photon
  Ctx c { t, par.nbin };
  staged_.resize(n);
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  thrust::fill(dbudget_.begin(), dbudget_.end(), 0.0);
  const auto stg  = thrust::raw_pointer_cast(staged_.data());
  const auto cnt  = thrust::raw_pointer_cast(counters_.data());
  const auto db   = thrust::raw_pointer_cast(dbudget_.data());
  const auto key  = runko::rng_key(par.seed, istep_, 0, 303);
  const auto mins = t.mins, maxs = t.maxs;                       // tile extent for the positions
  const auto lay  = c.v[0].lay;
  const auto w    = tyvi::mdgrid_work {};
  const auto i0   = thrust::counting_iterator<I>(0);
  // one lane per photon: a uniform position, an energy from the spectrum, an isotropic direction
  thrust::for_each(w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
    auto rng = runko::Rng { key, i };
    const float pos[3] = { mins[0] + (maxs[0] - mins[0]) * static_cast<float>(rng.uniform()),
                           mins[1] + (maxs[1] - mins[1]) * static_cast<float>(rng.uniform()),
                           mins[2] + (maxs[2] - mins[2]) * static_cast<float>(rng.uniform()) };
    // Wien: a Gamma(3, Theta) deviate; Planck: the term j of the Bose series first, then Gamma(3)/j
    auto j = 1.0;
    if(planck) {
      const auto z1 = 1.2020569 * rng.uniform();
      auto cum      = 1.0;
      while(z1 > cum and j < 50.0) {
        j += 1.0;
        cum += 1.0 / (j * j * j);
      }
    }
    const auto x    = -theta * sstd::log(rng.uniform() * rng.uniform() * rng.uniform()) / j;
    const auto k    = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    stg[k]          = staged_at(x * iso_dir(rng), pos, wph, 1, lay.cell_of(pos[0], pos[1], pos[2]), 0, static_cast<int>(rng.uniform() < 0.5));
    sstd::atomic_add(&db[b_E], wph * x);
  });
  w.wait();
  {
    // the energy injected
    thrust::host_vector<double> h = dbudget_;
    budget.inj_E += h[b_E];
  }
  this->place(c);
}

void
  Sampler::heat_diffusive(TileView& t, const double ell)
{
  this->lay_if_needed(t);
  // The Ito kick of the isotropic Fokker-Planck operator (1/p^2) d/dp [p^2 D dn/dp], taken on the
  // momentum vector so that it also scatters the direction. Its expected energy gain per lepton is
  // D_0 p^2 (4 gamma^2 + 1)/gamma^3, and D_0 is chosen every step so that the summed gain is exactly
  // the energy the compactness injects.
  Ctx c { t, par.nbin };
  const auto w  = tyvi::mdgrid_work {};
  auto tot      = 0.0;   // summed expected gain per unit D_0 over both lepton species
  // ---- the normalization: reduce w p^2 (4 gamma^2 + 1)/gamma^3 over the live leptons
  for(auto s = 1uz; s <= 2uz; ++s) {
    const auto v  = c.v[s];
    const auto n  = t.con[s]->size();
    const auto i0 = thrust::counting_iterator<I>(0);
    tot += thrust::transform_reduce(
      w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return 0.0; }
        const auto p   = v.u(i);
        const auto p2  = toolbox::dot(p, p);
        const auto gam = sstd::sqrt(1.0 + p2);
        return v.w(i) * p2 * (4.0 * gam * gam + 1.0) / (gam * gam * gam);
      },
      0.0, thrust::plus<double> {});
  }
  if(tot <= 0.0) { return; }
  const auto D0 = box.dE(ell) / tot;   // the diffusion coefficient that injects exactly dE
  thrust::fill(dbudget_.begin(), dbudget_.end(), 0.0);
  const auto db = thrust::raw_pointer_cast(dbudget_.data());   // energy gained
  // ---- the kicks: one lane per lepton slot of either species
  for(auto s = 1uz; s <= 2uz; ++s) {
    const auto v   = c.v[s];
    const auto n   = t.con[s]->size();
    const auto key = runko::rng_key(par.seed, istep_, s, 304);
    const auto i0  = thrust::counting_iterator<I>(0);
    thrust::for_each(w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      if(not v.alive(i)) { return; }
      auto rng      = runko::Rng { key, i };
      const auto p  = v.u(i);
      const auto p2 = toolbox::dot(p, p);
      const auto sd = sstd::sqrt(2.0 * D0 * p2);        // width of the Ito kick
      const auto p1 = (1.0 + 2.0 * D0) * p + Vec3(sd * rng.normal(), sd * rng.normal(), sd * rng.normal());   // drift plus noise
      sstd::atomic_add(&db[b_E], v.w(i) * (sstd::sqrt(1.0 + toolbox::dot(p1, p1)) - sstd::sqrt(1.0 + p2)));
      v.set_u(i, p1);
      v.touch(i);
    });
  }
  w.wait();
  // the energy the heating put in
  thrust::host_vector<double> h = dbudget_;
  budget.heat_E += h[b_E];
}

void
  Sampler::escape_photons(TileView& t, const double tau_ext, const double dt_ext)
{
  this->lay_if_needed(t);
  Ctx c { t, par.nbin };
  // the pair depth the escape sees: the whole slab's, given by the driver, or this tile's own
  const auto tau_T  = tau_ext >= 0.0 ? tau_ext : Pop::totals(*t.con[1]).W + Pop::totals(*t.con[2]).W;
  const auto dt_esc = dt_ext > 0.0 ? dt_ext : box.dt;   // the step in crossings of that slab
  if(esc_hist_.size() != static_cast<std::size_t>(NESC)) {
    esc_hist_.assign(NESC, 0.0);
    esc_mu_hist_.assign(static_cast<std::size_t>(NESC) * NMU, 0.0);
  }
  thrust::fill(dbudget_.begin(), dbudget_.end(), 0.0);
  const auto db   = thrust::raw_pointer_cast(dbudget_.data());   // escaped energy and weight
  const auto hist = thrust::raw_pointer_cast(esc_hist_.data());  // escaping x ell_x per log bin
  const auto hmu  = thrust::raw_pointer_cast(esc_mu_hist_.data());   // and per (log bin, direction)
  const auto v    = c.v[0];
  const auto n    = t.con[0]->size();
  const auto key  = runko::rng_key(par.seed, istep_, 0, 305);
  const auto llo  = std::log(ESC_LO), dl = (std::log(ESC_HI) - std::log(ESC_LO)) / NESC;   // the log bins
  const auto w    = tyvi::mdgrid_work {};
  const auto i0   = thrust::counting_iterator<I>(0);
  // one lane per photon slot: escape with the slab probability, book the energy
  thrust::for_each(w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
    if(not v.alive(i)) { return; }
    auto rng     = runko::Rng { key, i };
    const auto x = v.e(i);
    if(not(rng.uniform() < Box::p_esc(x, tau_T, dt_esc))) { return; }
    const auto wg = v.w(i);
    v.kill(i);
    sstd::atomic_add(&db[b_E], x * wg);
    sstd::atomic_add(&db[b_W], wg);
    const auto b = sstd::clamp(static_cast<int>(sstd::floor((sstd::log(x) - llo) / dl)), 0, NESC - 1);   // its spectral bin
    sstd::atomic_add(&hist[b], x * wg);
    const auto mu = static_cast<double>(v.vel[i][2]) / x;   // escape direction, cosine to z
    const auto m  = sstd::clamp(static_cast<int>(sstd::floor((mu + 1.0) * 0.5 * NMU)), 0, NMU - 1);
    sstd::atomic_add(&hmu[b * NMU + m], x * wg);
  });
  w.wait();
  // the escaped energy and weight
  thrust::host_vector<double> h = dbudget_;
  budget.esc_E += h[b_E];
  budget.esc_W += h[b_W];
}

// host copy of the escaping spectrum per direction, zeros before the first escape
std::vector<double>
  Sampler::escaped_mu_spectrum() const
{
  if(esc_mu_hist_.empty()) { return std::vector<double>(static_cast<std::size_t>(NESC) * NMU, 0.0); }
  thrust::host_vector<double> h = esc_mu_hist_;
  return std::vector<double>(h.begin(), h.end());
}

// host copy of the escaping spectrum, zeros before the first escape
std::vector<double>
  Sampler::escaped_spectrum() const
{
  if(esc_hist_.empty()) { return std::vector<double>(NESC, 0.0); }
  thrust::host_vector<double> h = esc_hist_;
  return std::vector<double>(h.begin(), h.end());
}

// start a new accounting window: budgets and the escaping spectrum
void
  Sampler::reset_budget()
{
  budget = Budget {};
  thrust::fill(esc_hist_.begin(), esc_hist_.end(), 0.0);
  thrust::fill(esc_mu_hist_.begin(), esc_mu_hist_.end(), 0.0);
}

double
  Sampler::tau_ssa(const TileView& t, const double x) const
{
  // With weight in units of depth, the absorptions per crossing of a photon are the sum over the
  // leptons of their weight times 2 sigma-hat F: the depth of the box in this channel, the analogue
  // of tau_T. Its unit crossing marks x_SSA.
  if(box.b_perp() <= 0.0) { return 0.0; }
  const auto fv = t.field;             // each lepton absorbs at the tile field across its pitch angle
  const auto w  = tyvi::mdgrid_work {};
  auto out      = 0.0;
  const auto k  = Vec3(x, 0.0, 0.0);   // the rate depends on x, gamma and the field only
  // sum w 2 sigma-hat F over the live leptons of both species
  for(auto s = 1uz; s <= 2uz; ++s) {
    auto& con     = *t.con[s];
    const auto v  = make_view(con);
    const auto n  = con.size();
    const auto i0 = thrust::counting_iterator<I>(0);
    out += thrust::transform_reduce(
      w.on_this(), i0, i0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return 0.0; }
        const auto u   = v.u(i);
        const auto b   = fv.b_ssa(u / sstd::max(toolbox::norm(u), 1e-300));
        return v.w(i) * 2.0 * ssa::rate(k, u, b, 0.0);
      },
      0.0, thrust::plus<double> {});
  }
  return out;
}

// ---------------------------------------------------------------- diagnostics
double
  Sampler::pair_sum(const Channel which, const TileView& t, const bool weight) const
{
  const auto w  = tyvi::mdgrid_work {};
  const auto N2 = box.N2();
  const auto fv = t.field;   // the partner's field, for the absorption rate
  // the brute-force sum over every (incident, partner) pair of the channel's record type
  auto run      = [&](const auto& ch) {
    const auto i1  = static_cast<std::size_t>(std::remove_cvref_t<decltype(ch)>::t1);
    const auto i2  = static_cast<std::size_t>(ch.t2);
    auto& c1       = *t.con[i1];
    auto& c2       = *t.con[i2];
    const auto v1 = make_view(c1);
    const auto v2 = make_view(c2);
    const auto n1 = c1.size(), n2 = c2.size();
    const auto r0 = thrust::counting_iterator<I>(0);
    // one lane per incident slot, looping over every partner slot
    auto tot      = thrust::transform_reduce(
      w.on_this(), r0, r0 + static_cast<std::ptrdiff_t>(n1),
      [=](const I i) {
        if(not v1.alive(i)) { return 0.0; }
        const auto p1 = v1.u(i);
        const auto wi = v1.w(i);
        auto s        = 0.0;
        for(I j = 0; j < n2; ++j) {
          if(not v2.alive(j)) { continue; }
          const auto wj = v2.w(j);
          s += sstd::max(wi, wj) * (weight ? sstd::min(wi, wj) : 1.0) * 2.0 * pair_rate(ch, fv, p1, v2.u(j));
        }
        return s;
      },
      0.0, thrust::plus<double> {});
    if(i1 == i2) { tot *= 0.5; }  // each unordered pair once
    return tot / N2;
  };
  // dispatch on the channel's record type
  switch(which) {
    case Channel::compton_m: return run(ChCompton { Species::em });
    case Channel::compton_p: return run(ChCompton { Species::ep });
    case Channel::pairann: return run(ChAnn {});
    case Channel::photann: return run(ChGG {});
    case Channel::ssa_m: return run(ChSSA { Species::em, fv.b_tile, box.x_min });
    case Channel::ssa_p: return run(ChSSA { Species::ep, fv.b_tile, box.x_min });
    default: throw std::runtime_error { "pair_sum: not a two-body channel" };
  }
}

// accepted events per enabled channel, cumulative
std::map<std::string, double>
  Sampler::counts() const
{
  std::map<std::string, double> m;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { m[NAMES[i]] = stats.ch[i].acc; }
  }
  return m;
}

// weight moved per enabled channel, cumulative: sum f min(w_i, w_j) over the events
std::map<std::string, double>
  Sampler::moved() const
{
  std::map<std::string, double> m;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { m[NAMES[i]] = stats.ch[i].wmoved; }
  }
  return m;
}

std::string
  Sampler::report() const
{
  std::ostringstream o;
  const auto ns = static_cast<double>(std::max<std::size_t>(1, stats.steps));   // everything per step
  o << std::format("kernels per step over {} steps ({} bin(s), sort every {}, slack {}):\n", stats.steps, par.nbin, par.nsort, par.slack);
  o << std::format("   slots / live LPs       : ph={:.0f}/{:.0f}, e-={:.0f}/{:.0f}, e+={:.0f}/{:.0f}\n", stats.slots[0] / ns, stats.live[0] / ns, stats.slots[1] / ns, stats.live[1] / ns, stats.slots[2] / ns, stats.live[2] / ns);
  o << std::format("   cell-kernel lanes      : {:.0f}\n", stats.cell_lanes / ns);
  // one line per enabled channel
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(not enabled[i]) { continue; }
    const auto& s = stats.ch[i];
    if(is_twobody(static_cast<Channel>(i))) {
      o << std::format("   {:9s}: lanes {:.1f}, virtual {:.1f}%, weight-stage {:.3f}, accepted {:.2f}, <sigmahat> {:.4f}, <f> {:.2f}, nsub {} ({} steps), clamped {:.0f}, trials {:.2f}\n",
        NAMES[i], s.M / ns, 100.0 * s.virt / std::max(1.0, s.M), s.wpass / std::max(1.0, s.M - s.virt), s.acc / ns,
        s.shbar_num / std::max(1e-30, s.shbar_den), s.facc / std::max(1.0, s.acc), s.nsub_max_seen, s.nsub_steps,
        s.clamped, s.trials / std::max(1.0, s.trials_n));
    } else {
      o << std::format("   {:9s}: lanes {:.0f}, events {:.2f} (dead {:.1f}%), emitted weight/step {:.3g}, <f> {:.2f}\n",
        NAMES[i], s.M / ns, s.acc / ns, 100.0 * s.virt / std::max(1.0, s.M), s.wem / ns, s.facc / std::max(1.0, s.acc));
    }
  }
  o << std::format("   products slack / tail  : {:.2f}/{:.2f}\n", stats.slackprod / ns, stats.tailprod / ns);
  return o.str();
}

}  // namespace qed
