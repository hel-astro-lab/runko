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

constexpr const char* NAMES[] = { "compton-", "compton+", "pairann", "photann", "ssa-",
                                  "ssa+",     "syn-",     "syn+",    "bw" };

// ---------------------------------------------------------------- channel records
// The physics bound to records, nothing dispatched on names later: the sampler kernels are
// templates over these types.
struct ChCompton {
  Species t2;
  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::scatter;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  constexpr double rate(const Vec3& k, const Vec3& u) const { return comp::rate(k, u); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  {
    return comp::bound(lo1, hi1, lo2, hi2);
  }
  constexpr double dgmax(double x, double g) const { return comp::dgmax(x, g); }
};

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
  constexpr double dgmax(double, double) const { return 0.0; }
};

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
  constexpr double dgmax(double, double) const { return 0.0; }
};

struct ChSSA {
  Species t2;
  double b_perp, x_min;
  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::absorb;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  constexpr double rate(const Vec3& k, const Vec3& u) const
  {
    return ssa::rate(k, u, b_perp, x_min);
  }
  constexpr double bound(double lo1, double /*hi1*/, double lo2, double hi2) const
  {
    return ssa::bound(b_perp, x_min, lo1, lo2, hi2);
  }
  constexpr double dgmax(double x, double) const { return ssa::dgmax(x, 0.0); }
};

// ---------------------------------------------------------------- device views
using VecMds    = decltype(std::declval<pic::ParticleContainer&>().vel_mds());
using IdMds     = decltype(std::declval<pic::ParticleContainer&>().ids_mds());
using WgtMds    = decltype(std::declval<pic::ParticleContainer&>().wgt_mds());

/// one species as the kernels see it
struct PView {
  VecMds pos, vel;
  IdMds ids;
  WgtMds wgt;
  bool weighted;
  float w0;
  float mass;
  LayoutView lay;

  [[nodiscard]] constexpr bool alive(const I i) const { return ids[i][] != runko::dead_prtc_id; }
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
  /// widen the bounds of the key of a slot updated in place
  constexpr void touch(const I i) const
  {
    const auto k = lay.skey[i];
    if(k < 0) { return; }
    const auto ee = static_cast<float>(e(i));
    const auto ww = static_cast<float>(w(i));
    sstd::atomic_max(&lay.emax[k], ee);
    sstd::atomic_min(&lay.emin[k], ee);
    sstd::atomic_max(&lay.wmax[k], ww);
  }
};

PView
  make_view(pic::ParticleContainer& con, Pop& pop)
{
  pop.follow_container(con);
  return PView { .pos      = con.pos_mds(),
                 .vel      = con.vel_mds(),
                 .ids      = con.ids_mds(),
                 .wgt      = con.wgt_mds(),
                 .weighted = con.has_particle_weights(),
                 .w0       = static_cast<float>(con.weight()),
                 .mass     = static_cast<float>(con.args().mass),
                 .lay      = pop.view() };
}

/// one side of a channel: the one or two sub-ranges of every (cell, bin) as one virtual range
struct SideView {
  const I* base0;
  const I* cnt0;
  const I* base1;
  const I* cnt1;
  const float* wmax;
  const float* emin;
  const float* emax;
  [[nodiscard]] constexpr I N(const int cb) const { return cnt0[cb] + cnt1[cb]; }
  /// virtual index v in [0, N) of group cb -> slot
  [[nodiscard]] constexpr I pick(const int cb, const I v) const
  {
    const auto c0 = cnt0[cb];
    return v < c0 ? base0[cb] + v : base1[cb] + v - c0;
  }
};

/// offset inside a range of N slots cut into M strata of a slot drawn uniformly from stratum p
[[nodiscard]]
constexpr I
  stratum(const std::uint64_t p, const std::uint64_t M, const std::uint64_t N, runko::Rng& rng)
{
  const auto B0 = (p * N) / M;
  const auto B1 = ((p + 1) * N) / M;
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

/// resolution function: LPs per unit weight go as (x/x_ref)^a for photons, 1 for leptons
[[nodiscard]]
constexpr double
  f_res(const Species t, const double e, const double a, const double x_ref)
{
  return t == Species::ph ? sstd::pow(e / x_ref, a) : 1.0;
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

enum DStat : int { d_virt = 0, d_wpass, d_wmoved, d_facc, d_trials, d_trials_n, d_clamped, d_overlap, d_wem, d_count };
enum Counter : int { c_staged = 0, c_tail0, c_tail1, c_tail2, c_slack, c_tailn, c_id0, c_id1, c_id2, c_count };

}  // namespace

// ---------------------------------------------------------------- the per-step context
struct Sampler::Ctx {
  TileView* t;
  std::array<PView, 3> v;
  std::array<std::size_t, 3> old_size {};
  int ncell;
  int nb;

  void refresh(const Species s)
  {
    const auto i = static_cast<std::size_t>(s);
    v[i]         = make_view(*t->con[i], *t->pop[i]);
  }
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
  const auto field = b.b_perp() > 0.0;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    const auto c    = static_cast<Channel>(i);
    const auto need = static_cast<int>(c) >= static_cast<int>(Channel::ssa_m);
    auto on         = not need or field;
    if(not channels.empty()) {
      on = on and std::ranges::find(channels, std::string(NAMES[i])) != channels.end();
    }
    enabled[i] = on;
  }
  counters_.assign(c_count, 0);
  dstat_.assign(d_count, 0.0);
}

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
    const std::string& site,
    const int ncell,
    const std::size_t n,
    const std::int32_t* cell,
    const double* q,
    std::int32_t* out)
{
  if(n == 0) { return; }
  auto& car        = carry(site, ncell);
  const auto car_p = thrust::raw_pointer_cast(car.data());
  const auto w     = tyvi::mdgrid_work {};
  const auto sorted =
    thrust::is_sorted(w.on_this(), cell, cell + static_cast<std::ptrdiff_t>(n));
  ev_cs_.resize(n);
  const auto cs = thrust::raw_pointer_cast(ev_cs_.data());
  const auto ord0 = thrust::counting_iterator<I>(0);
  if(sorted) {
    thrust::inclusive_scan_by_key(
      w.on_this(), cell, cell + static_cast<std::ptrdiff_t>(n), q, cs);
    thrust::for_each(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      const auto base = car_p[cell[i]];
      const auto m    = sstd::floor(base + cs[i]) - sstd::floor(base + cs[i] - q[i]);
      out[i]          = static_cast<std::int32_t>(m);
      sstd::atomic_add(&car_p[cell[i]], q[i] - m);
    });
  } else {
    // the tail of a species is not cell-ordered: sort a permutation by cell first
    ev_perm_.resize(n);
    auto keys = thrust::device_vector<std::int32_t>(cell, cell + static_cast<std::ptrdiff_t>(n));
    thrust::copy(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), ev_perm_.begin());
    thrust::stable_sort_by_key(w.on_this(), keys.begin(), keys.end(), ev_perm_.begin());
    const auto perm  = thrust::raw_pointer_cast(ev_perm_.data());
    const auto keys_p = thrust::raw_pointer_cast(keys.data());
    const auto qs    = thrust::make_transform_iterator(ord0, [=](const I i) { return q[perm[i]]; });
    thrust::inclusive_scan_by_key(
      w.on_this(), keys_p, keys_p + static_cast<std::ptrdiff_t>(n), qs, cs);
    thrust::for_each(w.on_this(), ord0, ord0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      const auto j    = perm[i];
      const auto base = car_p[keys_p[i]];
      const auto m    = sstd::floor(base + cs[i]) - sstd::floor(base + cs[i] - q[j]);
      out[j]          = static_cast<std::int32_t>(m);
      sstd::atomic_add(&car_p[keys_p[i]], q[j] - m);
    });
  }
  w.wait();
}

// ---------------------------------------------------------------- the step
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

void
  Sampler::step(TileView& t)
{
  auto need = istep_ % static_cast<std::uint64_t>(par.nsort) == 0;
  for(auto s = 0uz; s < 3uz; ++s) { need = need or not t.pop[s]->laid(); }
  if(need) { this->layout(t); }

  Ctx c { .t = &t, .v = {}, .old_size = {}, .ncell = 0, .nb = par.nbin };
  c.refresh_all();
  c.ncell = c.v[0].lay.ncell;
  for(auto s = 0uz; s < 3uz; ++s) {
    stats.slots[s] += static_cast<double>(t.con[s]->size());
    stats.live[s] += Pop::totals(*t.con[s]).N;
    t.ids_used[s] = 0;
  }
  thrust::fill(counters_.begin(), counters_.end(), 0);

  // The channels run one after another in a random order, each reading the states the
  // previous one left; the splitting error is second order in the per-channel rates and the
  // random order cancels it in expectation.
  auto order = std::vector<int> {};
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { order.push_back(static_cast<int>(i)); }
  }
  auto rng = runko::Rng { runko::rng_key(par.seed, istep_, 0, 0), 0 };
  for(auto i = order.size(); i > 1; --i) {
    const auto j = static_cast<std::size_t>(rng.integer(i));
    std::swap(order[i - 1], order[j]);
  }
  const auto b = box.b_perp();
  for(const auto l: order) {
    switch(static_cast<Channel>(l)) {
      case Channel::compton_m: twobody(c, ChCompton { Species::em }, l); break;
      case Channel::compton_p: twobody(c, ChCompton { Species::ep }, l); break;
      case Channel::pairann: twobody(c, ChAnn {}, l); break;
      case Channel::photann: twobody(c, ChGG {}, l); break;
      case Channel::ssa_m: twobody(c, ChSSA { Species::em, b, box.x_min }, l); break;
      case Channel::ssa_p: twobody(c, ChSSA { Species::ep, b, box.x_min }, l); break;
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
  const auto i1    = static_cast<std::size_t>(Ch::t1);
  const auto i2    = static_cast<std::size_t>(ch.t2);
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;
  const auto npair = ncb * nb;
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
    const auto b0 = thrust::raw_pointer_cast(base0.data()), c0 = thrust::raw_pointer_cast(cnt0.data());
    const auto wm = thrust::raw_pointer_cast(wmax.data());
    const auto lo = thrust::raw_pointer_cast(emin.data()), hi = thrust::raw_pointer_cast(emax.data());
    const auto inf = std::numeric_limits<float>::infinity();
    thrust::for_each(
      w.on_this(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(ncb), [=](const int cb) {
        int k0 = -1, k1 = -1;
        switch(sel) {
          case Sel::cell: k0 = cb; break;
          case Sel::bit0: k0 = 2 * cb; break;
          case Sel::bit1: k0 = 2 * cb + 1; break;
          case Sel::both: k0 = 2 * cb; k1 = 2 * cb + 1; break;
        }
        b0[cb] = lay.b[k0];
        c0[cb] = lay.fill[k0];
        auto wmx = lay.wmax[k0], emn = lay.emin[k0], emx = lay.emax[k0];
        if(k1 >= 0) {
          b0[ncb + cb] = lay.b[k1];
          c0[ncb + cb] = lay.fill[k1];
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
  const auto Lam = thrust::raw_pointer_cast(Lam_.data()), sh = thrust::raw_pointer_cast(sh_.data());
  const auto wmp = thrust::raw_pointer_cast(wm_.data());
  const auto rowsum = thrust::raw_pointer_cast(rowsum_.data()), colsum = thrust::raw_pointer_cast(colsum_.data());
  const auto kappa = par.kappa;
  const auto N2c   = box.N2() / ncell;
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
        // extremal at the corners of the gamma range, so they bound it from below
        auto wmv = sstd::max(w1, w2);
        if constexpr(Ch::has_dgm) {
          if(kappa > 0.0) {
            const auto fk = sstd::min(fkin(ch, kappa, hi1, lo2), fkin(ch, kappa, hi1, hi2));
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
  thrust::for_each(
    w.on_this(), thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(ncb), [=](const int cb) {
      const auto cc = cb / nb, bb = cb % nb;
      auto sum      = 0.0;
      for(auto a = 0; a < nb; ++a) { sum += Lam[(cc * nb + a) * nb + bb]; }
      colsum[cb] = sum;
    });
  stats.cell_lanes += npair;
  const auto g0 = thrust::counting_iterator<int>(0);
  st.shbar_num += thrust::transform_reduce(
    w.on_this(), g0, g0 + npair, [=](const int g) { return Lam[g] * sh[g]; }, 0.0, thrust::plus<double> {});
  st.shbar_den += thrust::transform_reduce(
    w.on_this(), g0, g0 + npair, [=](const int g) { return Lam[g]; }, 0.0, thrust::plus<double> {});

  // Distinct strata need at most as many events per bin as the bin has slots, on both sides;
  // where the envelope exceeds that the channel is split into sub-steps of a larger N_2. The
  // sub-steps of the additive channels are capped: beyond the cap the strata of a bin are let
  // overlap, so that a few heavy LPs are proposed as often as their rate demands.
  auto nsub = thrust::transform_reduce(
    w.on_this(), g0, g0 + ncb,
    [=](const int ca) {
      const auto n = static_cast<double>(s1.N(ca));
      return rowsum[ca] > 0.0 and n >= 1.0 ? sstd::ceil(rowsum[ca] / sstd::max(n - 1.0, 1.0)) : 1.0;
    },
    1.0, thrust::maximum<double> {});
  if constexpr(Ch::kind == Kind::convert) {
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
  for(auto r = 0; r < ns; ++r) { pass(c, ch, l, r, ns); }
}

template<typename Ch>
void
  Sampler::pass(Ctx& c, const Ch& ch, const int l, const int r, const int nsub)
{
  const auto i1    = static_cast<std::size_t>(Ch::t1);
  const auto i2    = static_cast<std::size_t>(ch.t2);
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;
  const auto npair = ncb * nb;
  auto& st         = stats.ch[static_cast<std::size_t>(l)];
  const auto w     = tyvi::mdgrid_work {};
  const auto step  = istep_;
  const auto keyA  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(1 + 4 * l));
  const auto keyB  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(2 + 4 * l));
  const auto keyC  = runko::rng_key(par.seed, step, static_cast<std::uint64_t>(r), static_cast<std::uint64_t>(3 + 4 * l));
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ds = thrust::raw_pointer_cast(dstat_.data());

  const auto s1 = SideView { thrust::raw_pointer_cast(base1_.data()), thrust::raw_pointer_cast(cnt1_.data()),
                             thrust::raw_pointer_cast(base1_.data()) + ncb, thrust::raw_pointer_cast(cnt1_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax1_.data()), thrust::raw_pointer_cast(emin1_.data()),
                             thrust::raw_pointer_cast(emax1_.data()) };
  const auto s2 = SideView { thrust::raw_pointer_cast(base2_.data()), thrust::raw_pointer_cast(cnt2_.data()),
                             thrust::raw_pointer_cast(base2_.data()) + ncb, thrust::raw_pointer_cast(cnt2_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax2_.data()), thrust::raw_pointer_cast(emin2_.data()),
                             thrust::raw_pointer_cast(emax2_.data()) };
  const auto Lam = thrust::raw_pointer_cast(Lam_.data()), sh = thrust::raw_pointer_cast(sh_.data());
  const auto wmp = thrust::raw_pointer_cast(wm_.data());

  // ---- cell kernel: counts by one cumulative rounding per incident bin, clamps, shifts
  M_.resize(static_cast<std::size_t>(npair));
  off1_.resize(static_cast<std::size_t>(npair));
  off2_.resize(static_cast<std::size_t>(npair));
  kk_.resize(static_cast<std::size_t>(npair));
  M1_.resize(static_cast<std::size_t>(ncb));
  M2_.resize(static_cast<std::size_t>(ncb));
  sh1_.resize(static_cast<std::size_t>(ncb));
  sh2_.resize(static_cast<std::size_t>(ncb));
  const auto M = thrust::raw_pointer_cast(M_.data()), off1 = thrust::raw_pointer_cast(off1_.data());
  const auto off2 = thrust::raw_pointer_cast(off2_.data()), kk = thrust::raw_pointer_cast(kk_.data());
  const auto M1 = thrust::raw_pointer_cast(M1_.data()), M2 = thrust::raw_pointer_cast(M2_.data());
  const auto sh1 = thrust::raw_pointer_cast(sh1_.data()), sh2 = thrust::raw_pointer_cast(sh2_.data());
  const auto cb0 = thrust::counting_iterator<int>(0);
  const auto convert = Ch::kind == Kind::convert;
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int ca) {
    auto rng  = runko::Rng { keyA, static_cast<std::uint64_t>(ca) };
    const auto U = rng.uniform();
    auto C = 0.0, Fp = 0.0;
    for(auto bb = 0; bb < nb; ++bb) {
      const auto g = ca * nb + bb;
      C += Lam[g] / nsub;
      const auto F = sstd::floor(C + U);
      M[g]         = static_cast<std::int32_t>(F - Fp);
      Fp           = F;
    }
    if(convert) {
      // hold the events of the row below its slot count: a backstop the sub-cycling makes rare
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
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int cb) {
    const auto cc = cb / nb, bb = cb % nb;
    if(convert) {
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
    if(not convert) {
      sstd::atomic_add(&ds[d_overlap], static_cast<double>(sstd::max(sum - static_cast<int>(s1.N(ca)), 0)));
    }
  });
  Mstart_.resize(static_cast<std::size_t>(npair) + 1);
  thrust::exclusive_scan(w.on_this(), M_.begin(), M_.end(), Mstart_.begin(), I { 0 });
  const auto Mtot = static_cast<std::size_t>(Mstart_[static_cast<std::size_t>(npair) - 1]) +
                    static_cast<std::size_t>(M_[static_cast<std::size_t>(npair) - 1]);
  Mstart_[static_cast<std::size_t>(npair)] = static_cast<I>(Mtot);
  st.M += static_cast<double>(Mtot);
  if(Mtot == 0) {
    w.wait();
    return;
  }

  // ---- thinning kernel: one lane per envelope event
  lane_acc_.resize(Mtot);
  lane_i_.resize(Mtot);
  lane_j_.resize(Mtot);
  lane_f_.resize(Mtot);
  const auto Mstart = thrust::raw_pointer_cast(Mstart_.data());
  const auto lacc = thrust::raw_pointer_cast(lane_acc_.data());
  const auto li = thrust::raw_pointer_cast(lane_i_.data()), lj = thrust::raw_pointer_cast(lane_j_.data());
  const auto lf = thrust::raw_pointer_cast(lane_f_.data());
  const auto v1 = c.v[i1], v2 = c.v[i2];
  const auto kappa = par.kappa;
  const auto lane0 = thrust::counting_iterator<I>(0);
  thrust::for_each(w.on_this(), lane0, lane0 + static_cast<std::ptrdiff_t>(Mtot), [=](const I lane) {
    auto rng     = runko::Rng { keyB, lane };
    const auto g = segment_of(Mstart, npair, lane);
    const auto e = static_cast<std::uint64_t>(lane - Mstart[g]);
    const auto cc = g / (nb * nb), a = (g / nb) % nb, bb = g % nb;
    const auto cb1 = cc * nb + a, cb2 = cc * nb + bb;
    const auto Mg = static_cast<std::uint64_t>(M[g]);
    const auto n1 = s1.N(cb1), n2 = s2.N(cb2);
    const auto p1 = static_cast<std::uint64_t>(off1[g]) + e;
    const auto vv1 = (static_cast<I>(sh1[cb1]) + stratum(p1, static_cast<std::uint64_t>(M1[cb1]), n1, rng)) % n1;
    I vv2;
    if(convert) {
      const auto p2 = static_cast<std::uint64_t>(off2[g]) + (e + static_cast<std::uint64_t>(kk[g])) % Mg;
      vv2 = (static_cast<I>(sh2[cb2]) + stratum(p2, static_cast<std::uint64_t>(M2[cb2]), n2, rng)) % n2;
    } else {
      // the partner only accumulates what the events give it, so it may be the partner of
      // several lanes: drawn uniformly with replacement
      vv2 = static_cast<I>(rng.uniform() * static_cast<double>(n2)) % n2;
    }
    const auto i = s1.pick(cb1, vv1), j = s2.pick(cb2, vv2);
    li[lane] = i;
    lj[lane] = j;
    lf[lane] = 1.0;
    lacc[lane] = 0;
    const auto virt = not(v1.alive(i) and v2.alive(j));
    if(virt) {
      sstd::atomic_add(&ds[d_virt], 1.0);
      return;
    }
    // Two acceptance stages against one uniform: the weight factor first, at the largest value
    // the merged pair can have; then the effective weight of the merged pair and the cross
    // section, whose threshold is the first times a factor <= 1.
    const auto u  = rng.uniform();
    const auto wi = v1.w(i), wj = v2.w(j);
    const auto wm = wmp[g];
    const auto wf = sstd::min(sstd::max(wi, wj) / wm, 1.0);
    if(not(u < wf)) { return; }
    sstd::atomic_add(&ds[d_wpass], 1.0);
    const auto p1v = v1.u(i), p2v = v2.u(j);
    const auto rate = ch.rate(p1v, p2v);
    auto weff = sstd::max(wi, wj);
    if constexpr(Ch::has_dgm) {
      if(kappa > 0.0) {
        const auto x   = toolbox::norm(p1v);
        const auto gam = sstd::sqrt(1.0 + toolbox::dot(p2v, p2v));
        weff           = sstd::max(wj, wi / fkin(ch, kappa, x, gam));
        lf[lane]       = sstd::max(wi, wj) / weff;
      }
    }
    if(u < (weff / wm) * rate / sh[g]) { lacc[lane] = 1; }
  });
  acc_.resize(Mtot);
  const auto acc_end = thrust::copy_if(
    w.on_this(), lane0, lane0 + static_cast<std::ptrdiff_t>(Mtot), acc_.begin(), [=](const I lane) { return lacc[lane] != 0; });
  const auto na = static_cast<std::size_t>(acc_end - acc_.begin());
  st.acc += static_cast<double>(na);
  {
    thrust::host_vector<double> h = dstat_;
    st.virt += h[d_virt];
    st.wpass += h[d_wpass];
    st.clamped += h[d_clamped];
    st.overlap += h[d_overlap];
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
  const auto acc = thrust::raw_pointer_cast(acc_.data());
  const auto q1 = thrust::raw_pointer_cast(ev_q1_.data()), q2 = thrust::raw_pointer_cast(ev_q2_.data());
  const auto q3 = thrust::raw_pointer_cast(ev_q3_.data());
  const auto n1p = thrust::raw_pointer_cast(ev_n1_.data()), n2p = thrust::raw_pointer_cast(ev_n2_.data());
  const auto n3p = thrust::raw_pointer_cast(ev_n3_.data());
  const auto ecell = thrust::raw_pointer_cast(ev_cell_.data());
  const auto ep = thrust::raw_pointer_cast(ev_p_.data());
  const auto wl = box.w_lep;
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ev0 = thrust::counting_iterator<I>(0);
  thrust::for_each(w.on_this(), ev0, ev0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
    auto rng       = runko::Rng { keyC, e };
    const auto lane = acc[e];
    const auto i = li[lane], j = lj[lane];
    const auto f = lf[lane];
    const auto g = segment_of(Mstart, npair, lane);
    ecell[e]     = g / (nb * nb);
    const auto p1v = v1.u(i), p2v = v2.u(j);
    const auto wi = v1.w(i), wj = v2.w(j);
    const auto s  = sstd::min(wi, wj);  // the weight one event transfers
    sstd::atomic_add(&ds[d_wmoved], f * s);
    sstd::atomic_add(&ds[d_facc], f);
    auto ntry = 0;
    Vec3 pa, pb;
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
    for(auto k = 0uz; k < 3uz; ++k) {
      ep[6 * e + k]     = static_cast<float>(pa[k]);
      ep[6 * e + 3 + k] = static_cast<float>(pb[k]);
    }
    sstd::atomic_add(&ds[d_trials], static_cast<double>(ntry));
    sstd::atomic_add(&ds[d_trials_n], 1.0);
  });
  w.wait();
  {
    thrust::host_vector<double> h = dstat_;
    st.wmoved += h[d_wmoved];
    st.facc += h[d_facc];
    st.trials += h[d_trials];
    st.trials_n += h[d_trials_n];
  }

  // ---- the decisions, by the systematic rounding of every cell
  const auto sp1 = std::string(NAMES[static_cast<int>(Ch::t1)]), sp2 = std::string(NAMES[static_cast<int>(ch.t2)]);
  std::ignore = sp1;
  std::ignore = sp2;
  const auto tag1 = std::to_string(static_cast<int>(Ch::t1)), tag2 = std::to_string(static_cast<int>(ch.t2));
  if constexpr(Ch::kind == Kind::scatter) {
    select("sc" + tag1, ncell, na, ecell, q1, n1p);
    select("kk" + tag2, ncell, na, ecell, q2, n2p);
  } else if constexpr(Ch::kind == Kind::absorb) {
    select("rm" + tag1, ncell, na, ecell, q1, n1p);
    select("kk" + tag2, ncell, na, ecell, q2, n2p);
  } else {
    select("rm" + tag1, ncell, na, ecell, q1, n1p);
    select("rm" + tag2, ncell, na, ecell, q2, n2p);
    if constexpr(Ch::t2 == Species::ph) { select("lep", ncell, na, ecell, q3, n3p); }
  }

  // ---- react kernel, second half: parents updated in place, products staged
  staged_.resize(2 * na);
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  touched_.assign(c.t->con[i2]->size(), std::uint8_t { 0 });
  const auto touched = thrust::raw_pointer_cast(touched_.data());
  if constexpr(Ch::kind == Kind::absorb) { dgam_.assign(c.t->con[i2]->size(), 0.0f); }
  const auto dgam = thrust::raw_pointer_cast(dgam_.data());
  const auto a_res = par.a_res, x_ref = par.x_ref;
  const auto n_cap = par.n_cap;
  const auto stage = [=](const I e, const Vec3& p, const float* pos, const double w_each, const int n_r, const int cell, const int species, const int bit) {
    const auto k = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    stg[k]       = Staged { .u       = { static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]) },
                            .pos     = { pos[0], pos[1], pos[2] },
                            .w       = static_cast<float>(w_each),
                            .n_r     = n_r,
                            .cell    = cell,
                            .species = species,
                            .bit     = bit };
    std::ignore  = e;
  };
  thrust::for_each(w.on_this(), ev0, ev0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
    auto rng        = runko::Rng { keyC, e + (std::uint64_t { 1 } << 32) };
    const auto lane = acc[e];
    const auto i = li[lane], j = lj[lane];
    const auto cell = ecell[e];
    const auto pa   = Vec3(static_cast<double>(ep[6 * e]), static_cast<double>(ep[6 * e + 1]), static_cast<double>(ep[6 * e + 2]));
    const auto pb   = Vec3(static_cast<double>(ep[6 * e + 3]), static_cast<double>(ep[6 * e + 4]), static_cast<double>(ep[6 * e + 5]));
    const float posi[3] = { v1.pos[i][0], v1.pos[i][1], v1.pos[i][2] };
    if constexpr(Ch::kind == Kind::scatter) {
      if(n1p[e] > 0) {
        // the parent takes the new state and is realized as the copies the resolution function
        // asks for; the parent is the first copy
        const auto x0 = v1.e(i), x1 = toolbox::norm(pa);
        const auto n  = f_res(Species::ph, x1, a_res, x_ref) / f_res(Species::ph, x0, a_res, x_ref);
        const auto cp = copies(n, v1.w(i), rng, n_cap);
        if(cp.n_r >= 1) {
          v1.set_u(i, pa);
          v1.set_w(i, cp.w_each);
          v1.touch(i);
          if(cp.n_r > 1) { stage(e, pa, posi, cp.w_each, cp.n_r - 1, cell, 0, v1.lay.skey[i] >= 0 ? v1.lay.skey[i] % 2 : 0); }
        } else {
          v1.kill(i);  // the roulette merged it away
        }
      }
      if(n2p[e] > 0) {
        for(auto k = 0uz; k < 3uz; ++k) { sstd::atomic_add(&v2.vel[j][k], static_cast<float>(pb[k])); }
        touched[j] = 1;
      }
    } else if constexpr(Ch::kind == Kind::absorb) {
      if(n2p[e] > 0) {
        sstd::atomic_add(&dgam[j], static_cast<float>(pa[0]));
        touched[j] = 1;
      }
      if(n1p[e] > 0) { v1.kill(i); }
    } else if constexpr(Ch::t2 == Species::ep) {
      const auto s  = sstd::min(v1.w(i), v2.w(j));
      const auto x0 = toolbox::norm(pa), x1 = toolbox::norm(pb);
      const auto ca = copies(f_res(Species::ph, x0, a_res, x_ref), s, rng, n_cap);
      if(ca.n_r >= 1) { stage(e, pa, posi, ca.w_each, ca.n_r, cell, 0, static_cast<int>(rng.uniform() < 0.5)); }
      const auto cb = copies(f_res(Species::ph, x1, a_res, x_ref), s, rng, n_cap);
      if(cb.n_r >= 1) { stage(e, pb, posi, cb.w_each, cb.n_r, cell, 0, static_cast<int>(rng.uniform() < 0.5)); }
      if(n1p[e] > 0) { v1.kill(i); }
      if(n2p[e] > 0) { v2.kill(j); }
    } else {
      const auto nl = n3p[e];
      if(nl > 0) {
        stage(e, pa, posi, wl, nl, cell, 1, 0);
        stage(e, pb, posi, wl, nl, cell, 2, 0);
      }
      if(n1p[e] > 0) { v1.kill(i); }
      if(n2p[e] > 0) { v2.kill(j); }
    }
  });
  w.wait();
  if constexpr(Ch::kind == Kind::absorb) { apply_dgam(c, ch.t2); }
  if constexpr(Ch::kind == Kind::scatter) { touch_marked(c, ch.t2); }
  this->place(c, keyC);
}

// the lepton takes the summed energy of its absorptions along its direction of motion
void
  Sampler::apply_dgam(Ctx& c, const Species s)
{
  const auto is = static_cast<std::size_t>(s);
  const auto v  = c.v[is];
  const auto dgam = thrust::raw_pointer_cast(dgam_.data());
  const auto touched = thrust::raw_pointer_cast(touched_.data());
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

void
  Sampler::touch_marked(Ctx& c, const Species s)
{
  const auto is = static_cast<std::size_t>(s);
  const auto v  = c.v[is];
  const auto touched = thrust::raw_pointer_cast(touched_.data());
  tyvi::mdgrid_work {}
    .for_each_index(v.ids, [=](const auto idx) {
      const auto i = static_cast<I>(idx[0]);
      if(touched[i] != 0 and v.alive(i)) { v.touch(i); }
    })
    .wait();
}

// ---------------------------------------------------------------- products: slack of their key, overflow to the tail
void
  Sampler::place(Ctx& c, const std::uint64_t key)
{
  const auto nst = static_cast<std::size_t>(counters_[c_staged]);
  if(nst == 0) { return; }
  const auto w   = tyvi::mdgrid_work {};
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto r0  = thrust::counting_iterator<I>(0);
  // copies per species, then room for them at the end of every container
  std::array<std::size_t, 3> tot {};
  for(auto s = 0; s < 3; ++s) {
    tot[static_cast<std::size_t>(s)] = static_cast<std::size_t>(thrust::transform_reduce(
      w.on_this(), r0, r0 + static_cast<std::ptrdiff_t>(nst),
      [=](const I k) { return stg[k].species == s ? static_cast<long>(stg[k].n_r) : 0L; }, 0L, thrust::plus<long> {}));
  }
  std::array<I, 3> tail_base {};
  for(auto s = 0uz; s < 3uz; ++s) {
    tail_base[s] = static_cast<I>(c.t->con[s]->size());
    if(tot[s] > 0) {
      c.t->con[s]->grow(tot[s]);
      c.t->pop[s]->follow_container(*c.t->con[s]);
    }
  }
  c.refresh_all();
  thrust::fill(counters_.begin() + c_tail0, counters_.end(), 0);
  const auto cnt = thrust::raw_pointer_cast(counters_.data());
  const auto v   = c.v;
  const auto idb = c.t->id_base;
  thrust::for_each(w.on_this(), r0, r0 + static_cast<std::ptrdiff_t>(nst), [=](const I k) {
    const auto rec = stg[k];
    const auto s   = static_cast<std::size_t>(rec.species);
    const auto pv  = v[s];
    const auto lay = pv.lay;
    const auto e   = sstd::sqrt(static_cast<double>(pv.mass) * pv.mass + static_cast<double>(rec.u[0]) * rec.u[0] +
                              static_cast<double>(rec.u[1]) * rec.u[1] + static_cast<double>(rec.u[2]) * rec.u[2]);
    const auto kk  = lay.key_of(rec.cell, lay.bin_of(e), rec.bit);  // the product's own key
    for(auto n = 0; n < rec.n_r; ++n) {
      const auto f = sstd::atomic_fetch_add(&lay.fill[kk], I { 1 });
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
    thrust::host_vector<std::int32_t> h = counters_;
    for(auto s = 0uz; s < 3uz; ++s) {
      c.t->ids_used[s] += static_cast<std::size_t>(h[c_id0 + static_cast<int>(s)]);
      c.t->id_base[s] += static_cast<runko::prtc_id_type>(h[c_id0 + static_cast<int>(s)]);
    }
    stats.slackprod[0] += h[c_slack];
    stats.tailprod[0] += h[c_tailn];
  }
  std::ignore = key;
}

// ---------------------------------------------------------------- one-body channels
void
  Sampler::onebody(Ctx& c, const Channel which, const int l)
{
  const auto emit = which != Channel::bw;
  const auto is   = static_cast<std::size_t>(which == Channel::syn_m ? Species::em : which == Channel::syn_p ? Species::ep : Species::ph);
  auto& st        = stats.ch[static_cast<std::size_t>(l)];
  const auto w    = tyvi::mdgrid_work {};
  c.refresh_all();
  const auto v    = c.v[is];
  const auto n    = c.t->con[is]->size();
  if(n == 0) { return; }
  const auto b    = box.b_perp();
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
  const auto keyA = runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(100 + 2 * l));
  const auto keyC = runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(101 + 2 * l));
  st.M += static_cast<double>(n);
  thrust::fill(dstat_.begin(), dstat_.end(), 0.0);
  const auto ds = thrust::raw_pointer_cast(dstat_.data());
  staged_.resize(2 * n);
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);
  const auto a_res = par.a_res, x_ref = par.x_ref, acc_syn = par.acc_syn, kappa = par.kappa;
  const auto n_cap = par.n_cap;
  const auto wl    = box.w_lep;
  const auto stage = [=](const Vec3& p, const I i, const double w_each, const int n_r, const int species, const int bit) {
    const auto k = static_cast<std::size_t>(sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    stg[k]       = Staged { .u       = { static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]) },
                            .pos     = { v.pos[i][0], v.pos[i][1], v.pos[i][2] },
                            .w       = static_cast<float>(w_each),
                            .n_r     = n_r,
                            .cell    = v.lay.cell_of(v.pos[i][0], v.pos[i][1], v.pos[i][2]),
                            .species = species,
                            .bit     = bit };
  };
  const auto s0 = thrust::counting_iterator<I>(0);
  if(emit) {
    // Below one emission per step the event is resolved: it fires with the merged probability and
    // emits one photon of the merged weight. Above it the parent cools continuously every step and
    // emits one photon LP per step standing for all the emissions of the step; the resolution
    // function then thins the soft emission by the roulette of copies().
    thrust::for_each(w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      if(not v.alive(i)) {
        sstd::atomic_add(&ds[d_virt], 1.0);
        return;
      }
      auto rngA     = runko::Rng { keyA, i };
      auto rngC     = runko::Rng { keyC, i };
      const auto u  = v.u(i);
      const auto lam = syn::rate(u, b, x_min, false) / N1;
      if(lam <= 0.0) { return; }
      const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
      const auto xm  = syn::xmean(gam * b, gam, false, x_min);
      const auto z   = sstd::max(rngA.uniform(), 1e-300);
      const auto slow = lam <= 1.0;
      const auto f    = slow ? facc_syn(gam, xm, acc_syn, kappa) : 1.0;
      const auto fire = slow ? -sstd::log(z) < lam / f : true;
      auto dgu = 0.0;
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
      const auto mult = slow ? 1.0 : dgu / sstd::max(xm, 1e-300);
      const auto P    = syn::sample(u, b, rngC, false, x_min);
      const auto x    = toolbox::norm(P.k);
      // a resolved parent gives up exactly what its photon carries, clipped at the mass shell; an
      // unresolved one has already cooled, and its LP carries the merged loss of f steps
      const auto dg  = slow ? f * x : dgu;
      const auto dgc = sstd::min(dg, gam - 1.0);
      const auto wem = v.w(i) * f * mult * (dg > 0.0 ? dgc / sstd::max(dg, 1e-300) : 1.0);
      if(slow) {
        const auto g1 = sstd::max(gam - dgc, 1.0);
        v.set_u(i, u * (sstd::sqrt(sstd::max(g1 * g1 - 1.0, 0.0)) / sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
      }
      v.touch(i);
      sstd::atomic_add(&ds[d_wpass], 1.0);  // fired
      sstd::atomic_add(&ds[d_wem], wem);
      sstd::atomic_add(&ds[d_facc], f);
      const auto cp = copies(f_res(Species::ph, x, a_res, x_ref), wem, rngC, n_cap);
      if(cp.n_r >= 1) { stage(P.k, i, cp.w_each, cp.n_r, 0, static_cast<int>(rngC.uniform() < 0.5)); }
    });
    w.wait();
    thrust::host_vector<double> h = dstat_;
    st.virt += h[d_virt];
    st.acc += h[d_wpass];
    st.wem += h[d_wem];
    st.facc += h[d_facc];
  } else {
    // magnetic pair production: fire, then the pair at the lepton weight, rounded systematically
    ev_q3_.resize(n);
    ev_n3_.resize(n);
    ev_cell_.resize(n);
    lane_acc_.assign(n, 0);
    const auto lacc = thrust::raw_pointer_cast(lane_acc_.data());
    thrust::for_each(w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), [=](const I i) {
      if(not v.alive(i)) {
        sstd::atomic_add(&ds[d_virt], 1.0);
        return;
      }
      auto rngA      = runko::Rng { keyA, i };
      const auto k   = v.u(i);
      const auto lam = bw::rate(k, b) / N1;
      const auto z   = sstd::max(rngA.uniform(), 1e-300);
      if(lam > 0.0 and -sstd::log(z) < lam) { lacc[i] = 1; }
    });
    acc_.resize(n);
    const auto acc_end = thrust::copy_if(
      w.on_this(), s0, s0 + static_cast<std::ptrdiff_t>(n), acc_.begin(), [=](const I i) { return lacc[i] != 0; });
    const auto na = static_cast<std::size_t>(acc_end - acc_.begin());
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
    thrust::for_each(w.on_this(), e0, e0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
      const auto i = acc[e];
      q3[e]        = v.w(i) / wl;
      ecell[e]     = v.lay.cell_of(v.pos[i][0], v.pos[i][1], v.pos[i][2]);
    });
    w.wait();
    select("lep", c.ncell, na, ecell, q3, n3);
    thrust::for_each(w.on_this(), e0, e0 + static_cast<std::ptrdiff_t>(na), [=](const I e) {
      const auto i = acc[e];
      auto rngC    = runko::Rng { keyC, i };
      const auto P = bw::sample(v.u(i), b, rngC);
      if(n3[e] > 0) {
        stage(P.um, i, wl, n3[e], 1, 0);
        stage(P.up, i, wl, n3[e], 2, 0);
      }
      v.kill(i);  // s = w: the parent always goes
    });
    w.wait();
  }
  this->place(c, keyC);
}

// ---------------------------------------------------------------- diagnostics
double
  Sampler::pair_sum(const Channel which, const TileView& t, const bool weight) const
{
  const auto w  = tyvi::mdgrid_work {};
  const auto N2 = box.N2();
  auto run      = [&](const auto& ch) {
    const auto i1  = static_cast<std::size_t>(std::remove_cvref_t<decltype(ch)>::t1);
    const auto i2  = static_cast<std::size_t>(ch.t2);
    auto& c1       = *t.con[i1];
    auto& c2       = *t.con[i2];
    const auto v1 = PView { c1.pos_mds(), c1.vel_mds(), c1.ids_mds(), c1.wgt_mds(), c1.has_particle_weights(), static_cast<float>(c1.weight()), static_cast<float>(c1.args().mass), {} };
    const auto v2 = PView { c2.pos_mds(), c2.vel_mds(), c2.ids_mds(), c2.wgt_mds(), c2.has_particle_weights(), static_cast<float>(c2.weight()), static_cast<float>(c2.args().mass), {} };
    const auto n1 = c1.size(), n2 = c2.size();
    const auto r0 = thrust::counting_iterator<I>(0);
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
          s += sstd::max(wi, wj) * (weight ? sstd::min(wi, wj) : 1.0) * 2.0 * ch.rate(p1, v2.u(j));
        }
        return s;
      },
      0.0, thrust::plus<double> {});
    if(i1 == i2) { tot *= 0.5; }  // each unordered pair once
    return tot / N2;
  };
  const auto b = box.b_perp();
  switch(which) {
    case Channel::compton_m: return run(ChCompton { Species::em });
    case Channel::compton_p: return run(ChCompton { Species::ep });
    case Channel::pairann: return run(ChAnn {});
    case Channel::photann: return run(ChGG {});
    case Channel::ssa_m: return run(ChSSA { Species::em, b, box.x_min });
    case Channel::ssa_p: return run(ChSSA { Species::ep, b, box.x_min });
    default: throw std::runtime_error { "pair_sum: not a two-body channel" };
  }
}

std::map<std::string, double>
  Sampler::counts() const
{
  std::map<std::string, double> m;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { m[NAMES[i]] = stats.ch[i].acc; }
  }
  return m;
}

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
  const auto ns = static_cast<double>(std::max<std::size_t>(1, stats.steps));
  o << std::format("kernels per step over {} steps ({} bin(s), sort every {}, slack {}):\n", stats.steps, par.nbin, par.nsort, par.slack);
  o << std::format("   slots / live LPs       : ph={:.0f}/{:.0f}, e-={:.0f}/{:.0f}, e+={:.0f}/{:.0f}\n", stats.slots[0] / ns, stats.live[0] / ns, stats.slots[1] / ns, stats.live[1] / ns, stats.slots[2] / ns, stats.live[2] / ns);
  o << std::format("   cell-kernel lanes      : {:.0f}\n", stats.cell_lanes / ns);
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(not enabled[i]) { continue; }
    const auto& s = stats.ch[i];
    if(is_twobody(static_cast<Channel>(i))) {
      o << std::format("   {:9s}: lanes {:.1f}, virtual {:.1f}%, weight-stage {:.3f}, accepted {:.2f}, <sigmahat> {:.4f}, <f> {:.2f}, nsub {} ({} steps), clamped {:.0f}, overlap {:.1f}%, trials {:.2f}\n",
        NAMES[i], s.M / ns, 100.0 * s.virt / std::max(1.0, s.M), s.wpass / std::max(1.0, s.M - s.virt), s.acc / ns,
        s.shbar_num / std::max(1e-30, s.shbar_den), s.facc / std::max(1.0, s.acc), s.nsub_max_seen, s.nsub_steps,
        s.clamped, 100.0 * s.overlap / std::max(1.0, s.M), s.trials / std::max(1.0, s.trials_n));
    } else {
      o << std::format("   {:9s}: lanes {:.0f}, events {:.2f} (dead {:.1f}%), emitted weight/step {:.3g}, <f> {:.2f}\n",
        NAMES[i], s.M / ns, s.acc / ns, 100.0 * s.virt / std::max(1.0, s.M), s.wem / ns, s.facc / std::max(1.0, s.acc));
    }
  }
  o << std::format("   products slack / tail  : {:.2f}/{:.2f}\n", stats.slackprod[0] / ns, stats.tailprod[0] / ns);
  return o.str();
}

}  // namespace qed
