// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

// The two-body channels. One channel on one tile runs as
//
//   sides -> envelope -> sub-steps of { strata -> lanes -> decide -> react -> settle }
//
// sides     per (cell, bin) group of both species: its slots and its bounds
// envelope  per (cell, bin a, bin b): an upper bound of the summed event probability
// strata    envelope events per bin pair, numbered row by row into lanes
// lanes     every lane proposes one slot pair and thins it to the pair's probability
// decide    what an accepted event does to its parents, rounded per cell
// react     the final state: incident in place, products staged, partner read only
// settle    the partners take what the events of the pass gave them

#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/sampler.h"
#include "runko/qed/weights.h"
#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "thrust/copy.h"
#include "thrust/for_each.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/scan.h"
#include "thrust/transform_reduce.h"
#include "tyvi/mdgrid.h"

#include <concepts>
#include <limits>

namespace qed {

namespace {

/// Returns a slot drawn uniformly from stratum p of N slots cut into M strata.
[[nodiscard]]
constexpr I
  stratum(
    const std::uint64_t p,
    const std::uint64_t M,
    const std::uint64_t N,
    runko::Rng& rng)
{
  // stratum p is [B0, B1)
  const auto B0 = (p * N) / M;
  const auto B1 = ((p + 1) * N) / M;

  return static_cast<I>(
    B0 + static_cast<std::uint64_t>(rng.uniform() * static_cast<double>(B1 - B0)));
}

/// Returns the index of the last start <= v in the sorted start[0, n).
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

// outcome of a thinning lane: weight stage failed, accepted, virtual (a dead or unkeyed
// slot), weight stage passed; the counts come from these, not from contended atomics
enum LaneState : std::uint8_t { L_WFAIL = 0, L_ACC = 1, L_VIRT = 2, L_WPASS = 3 };

// ---------------------------------------------------------------- reactions

/// An accepted event as its react lane holds it.
struct [[nodiscard]] Event {
  PView v1, v2;             // incident and partner species
  I i, j;                   // their slots
  double f;                 // exchanges the partner takes as one
  std::int32_t n1, n2, n3;  // rounded decisions: incident, partner, pair count
  int cell;
  float pos[3];  // products are born where the incident is
  Stager stage;
  double a_res, x_ref;  // photon resolution function
  int n_cap;            // largest copy number
  double w_lep;
  float* kick;  // the partner's momentum kick of this event
  float* dgam;  // energy absorbed, per partner slot
};

/// Compton: the photon takes the scattered state, as the copies the resolution function
/// asks for; the lepton's f exchanges are kept as one momentum kick.
constexpr void
  reaction(const ChCompton&, const Event& r, runko::Rng& rng)
{
  const auto p1v = r.v1.u(r.i), p2v = r.v2.u(r.j);
  auto ntry      = 0;  // proposals of the rejection sampler
  const auto P   = comp::sample(p1v, p2v, rng, ntry);
  const auto dp  = (-r.f) * (P.k1 - p1v);
  for(auto k = 0uz; k < 3uz; ++k) { r.kick[k] = static_cast<float>(dp[k]); }

  if(r.n1 <= 0) { return; }

  // the parent is the first copy
  const auto x0 = r.v1.e(r.i), x1 = toolbox::norm(P.k1);
  const auto n  = f_res(x1, r.a_res, r.x_ref) / f_res(x0, r.a_res, r.x_ref);
  const auto cp = copies(n, r.v1.w(r.i), rng, r.n_cap);
  if(cp.n_r < 1) {
    r.v1.kill(r.i);  // the roulette merged it away
    return;
  }

  // read before the touch, which drops a parent that left its bin from its key
  const auto bit = r.v1.lay.skey[r.i] % 2;
  r.v1.set_u(r.i, P.k1);
  r.v1.set_w(r.i, cp.w_each);
  r.v1.touch(r.i);
  if(cp.n_r > 1) {
    r.stage(P.k1, r.pos, cp.w_each, cp.n_r - 1, r.cell, Species::ph, bit);
  }
}

/// Self-absorption: the photon goes, the lepton gains the energy of f photons.
constexpr void
  reaction(const ChSSA&, const Event& r, runko::Rng&)
{
  if(r.n2 > 0) {
    // summed per lepton: other lanes may feed the same one
    const auto x = toolbox::norm(r.v1.u(r.i));
    sstd::atomic_add(&r.dgam[r.j], static_cast<float>(r.f * x));
  }
  if(r.n1 > 0) { r.v1.kill(r.i); }
}

/// Pair annihilation: two photons of the weight moved, each as the copies its
/// resolution asks for.
constexpr void
  reaction(const ChAnn&, const Event& r, runko::Rng& rng)
{
  auto ntry    = 0;
  const auto P = ann::sample(r.v1.u(r.i), r.v2.u(r.j), rng, ntry);
  const auto s = sstd::min(r.v1.w(r.i), r.v2.w(r.j));

  const auto photon = [&](const Vec3& k) {
    const auto cp  = copies(f_res(toolbox::norm(k), r.a_res, r.x_ref), s, rng, r.n_cap);
    const auto bit = static_cast<int>(rng.uniform() < 0.5);
    if(cp.n_r >= 1) { r.stage(k, r.pos, cp.w_each, cp.n_r, r.cell, Species::ph, bit); }
  };
  photon(P.k0);
  photon(P.k1);

  if(r.n1 > 0) { r.v1.kill(r.i); }
  if(r.n2 > 0) { r.v2.kill(r.j); }
}

/// Photon-photon pair production: the weight moved becomes n3 pairs at the lepton
/// weight.
constexpr void
  reaction(const ChGG&, const Event& r, runko::Rng& rng)
{
  auto ntry       = 0;
  const auto P    = gg::sample(r.v1.u(r.i), r.v2.u(r.j), rng, ntry);
  const auto flip = rng.uniform() < 0.5;

  if(r.n3 > 0) {
    // a pair born on top of itself adds no charge density: div E is unchanged
    r.stage(flip ? P.p1 : P.p2, r.pos, r.w_lep, r.n3, r.cell, Species::em, 0);
    r.stage(flip ? P.p2 : P.p1, r.pos, r.w_lep, r.n3, r.cell, Species::ep, 0);
  }
  if(r.n1 > 0) { r.v1.kill(r.i); }
  if(r.n2 > 0) { r.v2.kill(r.j); }
}

}  // namespace

// ---------------------------------------------------------------- one channel run

/// Channel ch (number l) on the tile of c: the stages of the header comment, over the
/// scratch of the sampler. Members end in an underscore so that a kernel cannot capture
/// one by accident: every stage copies what its lambdas read into locals.
template<typename Ch>
class TwoBody {
  Sampler& sam_;
  Ctx& c_;
  const Ch ch_;
  const Channel which_;
  const int l_;

  const int nb_;            // energy bins per cell
  const int ncell_;         // reaction cells
  const int ncb_;           // (cell, bin) groups
  const int npair_;         // (cell, bin a, bin b) pairs
  const std::size_t i1_;    // incident species
  const std::size_t i2_;    // partner species
  SideView s1_ {}, s2_ {};  // their sides, as of the start of the channel

  ChannelStats& st_;
  double* const ds_;  // the channel's device statistics

  static inline const auto g0   = thrust::counting_iterator<int>(0);  // groups, pairs
  static inline const auto ord0 = thrust::counting_iterator<I>(0);    // lanes, events

public:
  TwoBody(Sampler& sam, Ctx& c, const Ch& ch, const Channel which) :
    sam_ { sam },
    c_ { c },
    ch_ { ch },
    which_ { which },
    l_ { static_cast<int>(which) },
    nb_ { c.nb },
    ncell_ { c.ncell },
    ncb_ { c.ncell * c.nb },
    npair_ { c.ncell * c.nb * c.nb },
    i1_ { static_cast<std::size_t>(Ch::t1) },
    i2_ { static_cast<std::size_t>(ch.t2) },
    st_ { sam.stats.ch[static_cast<std::size_t>(which)] },
    ds_ { ptr(sam.dstat_) + static_cast<std::size_t>(which) * Sampler::d_count }
  {
  }

  void run()
  {
    const auto w = tyvi::mdgrid_work {};
    c_.refresh_all();

    s1_ = this->side(w, sam_.sides_[0], i1_, Ch::side1);
    s2_ = this->side(w, sam_.sides_[1], i2_, Ch::side2);
    this->envelope(w);
    const auto nsub = this->substeps(w);

    // each sub-step with the envelope divided by nsub
    for(auto r = 0; r < nsub; ++r) { this->pass(r, nsub); }
  }

private:
  /// Fills table t with the groups of species is and returns its view: the slot
  /// sub-ranges of every (cell, bin) and the bounds over them.
  [[nodiscard]]
  SideView side(
    const tyvi::mdgrid_work& w,
    Sampler::SideTable& t,
    const std::size_t is,
    const Side keys)
  {
    const auto ncb = ncb_;
    t.base.resize(static_cast<std::size_t>(2 * ncb));
    t.cnt.resize(static_cast<std::size_t>(2 * ncb));
    t.wmax.resize(static_cast<std::size_t>(ncb));
    t.emin.resize(static_cast<std::size_t>(ncb));
    t.emax.resize(static_cast<std::size_t>(ncb));

    const auto lay = c_.v[is].lay;
    const auto b0 = ptr(t.base), c0 = ptr(t.cnt);
    const auto wm = ptr(t.wmax), lo = ptr(t.emin), hi = ptr(t.emax);
    const auto inf = std::numeric_limits<float>::infinity();

    thrust::for_each(w.on_this(), g0, g0 + ncb, [=](const int cb) {
      // the group's one or two keys
      const auto k0 = keys.stride * cb + keys.first;
      const auto k1 = keys.nranges == 2 ? k0 + 1 : -1;

      // fill counts demand and exceeds cap after an overflow: only cap slots are
      // proposable, past them lies the next key or the tail
      b0[cb]   = lay.b[k0];
      c0[cb]   = sstd::min(lay.fill[k0], lay.cap[k0]);
      auto wmx = lay.wmax[k0], emn = lay.emin[k0], emx = lay.emax[k0];
      if(k1 >= 0) {
        b0[ncb + cb] = lay.b[k1];
        c0[ncb + cb] = sstd::min(lay.fill[k1], lay.cap[k1]);
        wmx          = sstd::max(wmx, lay.wmax[k1]);
        emn          = sstd::min(emn, lay.emin[k1]);
        emx          = sstd::max(emx, lay.emax[k1]);
      } else {
        b0[ncb + cb] = 0;
        c0[ncb + cb] = 0;
      }

      wm[cb] = wmx;
      lo[cb] = c0[cb] + c0[ncb + cb] > 0 ? emn : inf;
      hi[cb] = emx;
    });

    return SideView { b0, c0, b0 + ncb, c0 + ncb, wm, lo, hi };
  }

  /// Envelope per bin pair, from the slot counts and the bounds of the two bins alone:
  /// all n1 n2 slot pairs at the rate bound and the pair-weight bound of the bin pair.
  void envelope(const tyvi::mdgrid_work& w)
  {
    const auto nb = nb_, ncb = ncb_;
    auto& e = sam_.env_;
    e.Lam.resize(static_cast<std::size_t>(npair_));
    e.bound.resize(static_cast<std::size_t>(npair_));
    e.wmax.resize(static_cast<std::size_t>(npair_));
    e.row.resize(static_cast<std::size_t>(ncb));
    e.col.resize(static_cast<std::size_t>(ncb));

    const auto Lam = ptr(e.Lam), sh = ptr(e.bound), wmp = ptr(e.wmax);
    const auto rowsum = ptr(e.row), colsum = ptr(e.col);
    const auto s1 = s1_, s2 = s2_;
    const auto ch = ch_;
    const auto ds = ds_;

    const auto kappa = sam_.par.kappa;
    const auto N2c   = sam_.box.N2() / ncell_;  // the two-body clock of a reaction cell

    thrust::for_each(w.on_this(), g0, g0 + ncb, [=](const int ca) {
      const auto cc  = ca / nb;
      const auto n1  = static_cast<double>(s1.N(ca));
      const auto lo1 = static_cast<double>(s1.emin[ca]),
                 hi1 = static_cast<double>(s1.emax[ca]);
      const auto w1  = static_cast<double>(s1.wmax[ca]);
      auto sum = 0.0, sum_sh = 0.0;

      for(auto bb = 0; bb < nb; ++bb) {
        const auto cb  = cc * nb + bb;
        const auto g   = ca * nb + bb;
        const auto n2  = static_cast<double>(s2.N(cb));
        const auto lo2 = static_cast<double>(s2.emin[cb]),
                   hi2 = static_cast<double>(s2.emax[cb]);
        const auto w2  = static_cast<double>(s2.wmax[cb]);
        auto shv       = n1 > 0.0 and n2 > 0.0 ? ch.bound(lo1, hi1, lo2, hi2) : 0.0;

        // pair weight of the two maxima at the smallest f_kin of the bin pair, which
        // is extremal at the corners of the gamma range, starting at the threshold
        auto fk = 1.0;
        if constexpr(Ch::merges) {
          const auto glo = sstd::max(lo2, Ch::g_min);
          fk             = sstd::min(
            fkin(kappa, glo, ch.dgmax(hi1, glo)),
            fkin(kappa, hi2, ch.dgmax(hi1, hi2)));
        }
        const auto wmv = pair_weight(w1, w2, fk);

        auto L = Ch::c_l * n1 * n2 * 2.0 * shv * wmv / N2c;
        if(not(L == L) or L > 1e300) { L = 0.0; }

        Lam[g] = L;
        sh[g]  = shv;
        wmp[g] = wmv;
        sum += L;
        sum_sh += L * shv;
      }

      rowsum[ca] = sum;

      // the envelope-weighted mean of the rate bound, a measure of the bins' tightness
      sstd::atomic_add(&ds[Sampler::d_shnum], sum_sh);
      sstd::atomic_add(&ds[Sampler::d_shden], sum);
    });

    thrust::for_each(w.on_this(), g0, g0 + ncb, [=](const int cb) {
      const auto cc = cb / nb, bb = cb % nb;
      auto sum = 0.0;
      for(auto a = 0; a < nb; ++a) { sum += Lam[(cc * nb + a) * nb + bb]; }
      colsum[cb] = sum;
    });

    sam_.stats.cell_lanes += npair_;
  }

  /// Sub-steps of the channel. Distinct strata need at most as many events per bin as
  /// slots; beyond that the channel splits into sub-steps of a larger N_2. A consumed
  /// partner needs them on its side too. A channel that only adds to its partner caps
  /// them: past the cap the row clamp drops the surplus (an under-count the cap
  /// controls), so no two lanes ever hold one incident.
  [[nodiscard]]
  int substeps(const tyvi::mdgrid_work& w)
  {
    // sub-steps that bring the events of every group of side s below its slots
    const auto needed = [&](const SideView s, const double* const sum) {
      return thrust::transform_reduce(
        w.on_this(),
        g0,
        g0 + ncb_,
        [=](const int cb) {
          const auto n = static_cast<double>(s.N(cb));
          return sum[cb] > 0.0 and n >= 1.0
                   ? sstd::ceil(sum[cb] / sstd::max(n - 1.0, 1.0))
                   : 1.0;
        },
        1.0,
        thrust::maximum<double> {});
    };

    auto nsub = needed(s1_, ptr(sam_.env_.row));
    if constexpr(Ch::consumes_partner) {
      nsub = sstd::max(nsub, needed(s2_, ptr(sam_.env_.col)));
    } else {
      nsub = sstd::min(nsub, static_cast<double>(sam_.par.nsub_max));
    }

    const auto ns     = static_cast<int>(nsub);
    st_.nsub_max_seen = sstd::max(st_.nsub_max_seen, ns);
    st_.nsub_steps += (ns > 1);
    return ns;
  }

  /// Sub-step r of nsub: its lanes, its events and their reactions.
  void pass(const int r, const int nsub)
  {
    const auto w = tyvi::mdgrid_work {};

    // random streams of this (channel, sub-step): A the counts (lane = incident group)
    // and the shifts (lanes past them), B the thinning, C the reactions
    const auto key = [&](const int k) {
      return runko::rng_key(
        sam_.par.seed,
        sam_.istep_,
        static_cast<std::uint64_t>(r),
        stream::twobody(which_, k));
    };
    const auto keyA = key(0), keyB = key(1), keyC = key(2);

    const auto Mtot = this->strata(w, keyA, nsub);
    st_.M += static_cast<double>(Mtot);
    if(Mtot == 0) { return; }

    const auto na = this->lanes(w, keyA, keyB, Mtot);
    if(na == 0) {
      w.wait();
      return;
    }

    this->decide(w, na);
    this->react(w, keyC, na);
    w.wait();

    this->settle(na);
    sam_.place(c_);
  }

  /// Envelope events per bin pair for one sub-step, and the first lane of each; returns
  /// the lanes. The lanes run row by row, so a lane's rank in its row is its incident
  /// stratum; a consumed partner gets strata of its own along the columns.
  [[nodiscard]]
  std::size_t strata(const tyvi::mdgrid_work& w, const std::uint64_t keyA, const int nsub)
  {
    const auto nb = nb_, ncb = ncb_, npair = npair_;
    auto& t = sam_.strata_;
    t.M.resize(static_cast<std::size_t>(npair));
    t.off2.resize(static_cast<std::size_t>(npair));
    t.M2.resize(static_cast<std::size_t>(ncb));

    const auto M = ptr(t.M), off2 = ptr(t.off2), M2 = ptr(t.M2);
    const auto Lam = ptr(sam_.env_.Lam);
    const auto s1 = s1_, s2 = s2_;
    const auto ds = ds_;

    // incident rows: M[a, b] by cumulative rounding with one uniform for the whole row
    thrust::for_each(w.on_this(), g0, g0 + ncb, [=](const int ca) {
      auto rng     = runko::Rng { keyA, static_cast<std::uint64_t>(ca) };
      const auto U = rng.uniform();
      auto C = 0.0, Fp = 0.0;  // cumulative envelope and its previous floor

      for(auto bb = 0; bb < nb; ++bb) {
        const auto g = ca * nb + bb;
        C += Lam[g] / nsub;
        const auto F = sstd::floor(C + U);
        M[g]         = static_cast<std::int32_t>(F - Fp);
        Fp           = F;
      }

      // row clamp: no more events than slots, so the strata stay distinct; a backstop
      // the sub-steps make rare, counted as clamped
      auto ex = 0;
      for(auto bb = 0; bb < nb; ++bb) { ex += M[ca * nb + bb]; }
      ex -= static_cast<int>(s1.N(ca));

      for(; ex > 0; --ex) {
        auto imax = 0;
        for(auto bb = 1; bb < nb; ++bb) {
          if(M[ca * nb + bb] > M[ca * nb + imax]) { imax = bb; }
        }
        --M[ca * nb + imax];
        sstd::atomic_add(&ds[Sampler::d_clamped], 1.0);
      }
    });

    // partner columns of a consumed partner: the clamp, then the strata bases along a
    if constexpr(Ch::consumes_partner) {
      thrust::for_each(w.on_this(), g0, g0 + ncb, [=](const int cb) {
        const auto cc = cb / nb, bb = cb % nb;
        auto ex = 0;
        for(auto a = 0; a < nb; ++a) { ex += M[(cc * nb + a) * nb + bb]; }
        ex -= static_cast<int>(s2.N(cb));

        for(; ex > 0; --ex) {
          auto imax = 0;
          for(auto a = 1; a < nb; ++a) {
            if(M[(cc * nb + a) * nb + bb] > M[(cc * nb + imax) * nb + bb]) { imax = a; }
          }
          --M[(cc * nb + imax) * nb + bb];
          sstd::atomic_add(&ds[Sampler::d_clamped], 1.0);
        }

        // strata bases: disjoint across bin pairs
        auto sum = 0;
        for(auto a = 0; a < nb; ++a) {
          const auto g = (cc * nb + a) * nb + bb;
          off2[g]      = sum;
          sum += M[g];
        }
        M2[cb] = sum;
      });
    }

    // first lane of every bin pair; the last entry, the lane total, is read on the host
    t.first.resize(static_cast<std::size_t>(npair) + 1);
    t.first[0] = I { 0 };
    thrust::inclusive_scan(w.on_this(), t.M.begin(), t.M.end(), t.first.begin() + 1);
    w.wait();

    return static_cast<std::size_t>(t.first[static_cast<std::size_t>(npair)]);
  }

  /// The thinning: one lane per envelope event proposes a slot pair of its bin pair
  /// and accepts it with the pair's probability over the envelope's. Returns the
  /// accepted lanes, listed in acc_.
  [[nodiscard]]
  std::size_t lanes(
    const tyvi::mdgrid_work& w,
    const std::uint64_t keyA,
    const std::uint64_t keyB,
    const std::size_t Mtot)
  {
    const auto nb = nb_, ncb = ncb_, npair = npair_;
    sam_.lanes_.resize(Mtot);

    const auto rec   = ptr(sam_.lanes_);
    const auto first = ptr(sam_.strata_.first);
    const auto off2 = ptr(sam_.strata_.off2), M2 = ptr(sam_.strata_.M2);
    const auto sh = ptr(sam_.env_.bound), wmp = ptr(sam_.env_.wmax);
    const auto s1 = s1_, s2 = s2_;
    const auto ch = ch_;

    const auto v1 = c_.v[i1_], v2 = c_.v[i2_];
    const auto fv    = c_.t->field;  // the field at the partner, for the absorption rate
    const auto kappa = sam_.par.kappa;
    const auto nlane = static_cast<std::ptrdiff_t>(Mtot);

    thrust::for_each(w.on_this(), ord0, ord0 + nlane, [=](const I lane) {
      auto rng      = runko::Rng { keyB, lane };
      const auto g  = segment_of(first, npair, lane);
      const auto cc = g / (nb * nb), a = (g / nb) % nb,
                 bb  = g % nb;  // cell, incident bin, partner bin
      const auto cb1 = cc * nb + a, cb2 = cc * nb + bb;
      const auto n1 = s1.N(cb1), n2 = s2.N(cb2);

      // Every range is shifted cyclically by a uniform draw of its own, the same for
      // all lanes of the pass: whatever its stratum, a lane then lands on every slot
      // of the range with the same probability, independently on the two sides.
      const auto shift = [=](const int side, const int cb, const I n) {
        auto r = runko::Rng { keyA, static_cast<std::uint64_t>((1 + side) * ncb + cb) };
        return static_cast<I>(r.integer(static_cast<std::uint64_t>(n)));
      };

      // incident: a slot drawn in stratum p1 of the row's M1
      const auto row = first[cb1 * nb];
      const auto p1  = static_cast<std::uint64_t>(lane - row);
      const auto M1  = static_cast<std::uint64_t>(first[(cb1 + 1) * nb] - row);
      const auto vv1 = (shift(1, cb1, n1) + stratum(p1, M1, n1, rng)) % n1;

      I vv2;
      if constexpr(Ch::consumes_partner) {
        // a consumed partner: a slot drawn in stratum p2 of the column's M2
        const auto p2 = static_cast<std::uint64_t>(off2[g]) +
                        static_cast<std::uint64_t>(lane - first[g]);
        const auto M2c = static_cast<std::uint64_t>(M2[cb2]);
        vv2            = (shift(0, cb2, n2) + stratum(p2, M2c, n2, rng)) % n2;
      } else {
        // the partner only accumulates what the events give it, so several lanes may
        // share it: drawn uniformly with replacement
        vv2 = static_cast<I>(rng.uniform() * static_cast<double>(n2)) % n2;
      }

      const auto i = s1.pick(cb1, vv1), j = s2.pick(cb2, vv2);
      auto& ln = rec[lane];
      ln       = { .i = i, .j = j, .S = 0.0, .g = g, .state = L_WFAIL };

      // a dead slot, or one a kick dropped from its key: nothing to react
      if(not(v1.alive(i) and v2.alive(j)) or v1.lay.skey[i] < 0 or v2.lay.skey[j] < 0) {
        ln.state = L_VIRT;
        return;
      }

      // two stages on one uniform: the pair weight at its largest, without merging,
      // from the weights alone; then the pair weight and the rate from the momenta, a
      // threshold the first times a factor <= 1
      const auto u  = rng.uniform();
      const auto wi = v1.w(i), wj = v2.w(j);
      const auto wm = wmp[g];
      const auto wf = sstd::min(pair_weight(wi, wj, 1.0) / wm, 1.0);
      if(not(u < wf)) { return; }
      ln.state = L_WPASS;

      const auto p1v = v1.u(i), p2v = v2.u(j);
      const auto rate = ch.rate(p1v, p2v, fv);
      auto fk         = 1.0;
      if constexpr(Ch::merges) {
        const auto x   = toolbox::norm(p1v);
        const auto gam = sstd::sqrt(1.0 + toolbox::dot(p2v, p2v));
        fk             = fkin(kappa, gam, ch.dgmax(x, gam));
      }

      ln.S = event_size(wi, wj, fk);
      if(u < (pair_weight(wi, wj, fk) / wm) * rate / sh[g]) { ln.state = L_ACC; }
    });

    sam_.acc_.resize(Mtot);
    const auto acc_end =
      thrust::copy_if(w.on_this(), ord0, ord0 + nlane, sam_.acc_.begin(), [=](const I lane) {
        return rec[lane].state == L_ACC;
      });
    const auto na = static_cast<std::size_t>(acc_end - sam_.acc_.begin());
    st_.acc += static_cast<double>(na);

    // the dead lanes and those past the weight stage only, in one reduction
    using V2         = toolbox::VecD<double, 2>;
    const auto tally = thrust::transform_reduce(
      w.on_this(),
      ord0,
      ord0 + nlane,
      [=](const I lane) {
        return V2(
          rec[lane].state == L_VIRT ? 1.0 : 0.0,
          rec[lane].state == L_WPASS ? 1.0 : 0.0);
      },
      V2 {},
      [](const V2& a, const V2& b) { return a + b; });

    st_.virt += tally[0];
    st_.wpass += static_cast<double>(na) + tally[1];
    return na;
  }

  /// What every accepted event does to its parents. It follows from their weights and
  /// the size of the event alone, so it is settled before the final state is drawn:
  /// the probabilities, then their systematic rounding in every cell.
  void decide(const tyvi::mdgrid_work& w, const std::size_t na)
  {
    const auto nb = nb_;
    auto& ev = sam_.ev_;
    ev.q1.resize(na);
    ev.q2.resize(na);
    ev.q3.resize(na);
    ev.n1.resize(na);
    ev.n2.resize(na);
    ev.n3.resize(na);
    ev.cell.resize(na);
    ev.f.resize(na);

    const auto q1 = ptr(ev.q1), q2 = ptr(ev.q2), q3 = ptr(ev.q3), ef = ptr(ev.f);
    const auto ecell = ptr(ev.cell);
    const auto acc   = ptr(sam_.acc_);
    const auto rec   = ptr(sam_.lanes_);
    const auto v1 = c_.v[i1_], v2 = c_.v[i2_];
    const auto wl    = sam_.box.w_lep;
    const auto nev   = static_cast<std::ptrdiff_t>(na);

    thrust::for_each(w.on_this(), ord0, ord0 + nev, [=](const I e) {
      const auto ln = rec[acc[e]];
      const auto wi = v1.w(ln.i), wj = v2.w(ln.j);
      const auto S  = ln.S;    // the weight the event moves, <= wi
      const auto n2 = S / wj;  // the exchanges the partner takes
      ecell[e]      = ln.g / (nb * nb);

      // the incident gives up S; the partner takes n2 exchanges, with that probability
      // while n2 < 1 and as one event beyond
      q1[e] = S / wi;
      q2[e] = sstd::min(n2, 1.0);
      ef[e] = sstd::max(n2, 1.0);
      q3[e] = Ch::makes_pairs ? S / wl : 0.0;
    });
    // no wait: the rounding and the reactions follow on the same queue

    sam_.select(w, true, l_, 0, ncell_, na, ecell, q1, ptr(ev.n1));
    sam_.select(w, true, l_, 1, ncell_, na, ecell, q2, ptr(ev.n2));
    if constexpr(Ch::makes_pairs) {
      sam_.select(w, true, l_, 2, ncell_, na, ecell, q3, ptr(ev.n3));
    }
  }

  /// The final state of every event: the incident is updated in place and the products
  /// are staged. Lanes may share a partner, so its state is only read here; what the
  /// events give it is kept (a kick per event, absorbed energy per slot) for settle().
  void react(const tyvi::mdgrid_work& w, const std::uint64_t keyC, const std::size_t na)
  {
    auto& ev         = sam_.ev_;
    const auto stage = sam_.stager(2 * na);  // at most two products per event

    if constexpr(std::same_as<Ch, ChCompton>) { ev.kick.resize(3 * na); }

    // dgam is zero between passes (settle clears what it reads): only growth needs
    // new zeros
    if constexpr(std::same_as<Ch, ChSSA>) {
      if(sam_.dgam_.size() < c_.t->con[i2_]->size()) {
        sam_.dgam_.resize(c_.t->con[i2_]->size(), 0.0f);
      }
    }

    const auto n1p = ptr(ev.n1), n2p = ptr(ev.n2), n3p = ptr(ev.n3);
    const auto ecell = ptr(ev.cell);
    const auto ef    = ptr(ev.f);
    const auto kick  = ptr(ev.kick);
    const auto dgam  = ptr(sam_.dgam_);
    const auto acc   = ptr(sam_.acc_);
    const auto rec   = ptr(sam_.lanes_);
    const auto v1 = c_.v[i1_], v2 = c_.v[i2_];
    const auto ch    = ch_;
    const auto par   = sam_.par;
    const auto wl    = sam_.box.w_lep;
    const auto nev   = static_cast<std::ptrdiff_t>(na);

    thrust::for_each(w.on_this(), ord0, ord0 + nev, [=](const I e) {
      auto rng      = runko::Rng { keyC, e };
      const auto ln = rec[acc[e]];

      reaction(
        ch,
        Event { .v1    = v1,
                .v2    = v2,
                .i     = ln.i,
                .j     = ln.j,
                .f     = ef[e],
                .n1    = n1p[e],
                .n2    = n2p[e],
                .n3    = n3p[e],
                .cell  = ecell[e],
                .pos   = { v1.pos[ln.i][0], v1.pos[ln.i][1], v1.pos[ln.i][2] },
                .stage = stage,
                .a_res = par.a_res,
                .x_ref = par.x_ref,
                .n_cap = par.n_cap,
                .w_lep = wl,
                .kick  = kick + 3 * e,
                .dgam  = dgam },
        rng);
    });
  }

  /// The partners take what the events of the pass gave them, each once, and the
  /// statistics of the pass are read.
  void settle(const std::size_t na)
  {
    const auto w     = tyvi::mdgrid_work {};
    const auto v     = c_.v[i2_];
    const auto acc   = ptr(sam_.acc_);
    const auto rec   = ptr(sam_.lanes_);
    const auto ef    = ptr(sam_.ev_.f);
    const auto nev   = static_cast<std::ptrdiff_t>(na);

    // weight moved and exchanges merged
    using V2        = toolbox::VecD<double, 2>;
    const auto sums = thrust::transform_reduce(
      w.on_this(),
      ord0,
      ord0 + nev,
      [=](const I e) { return V2(rec[acc[e]].S, ef[e]); },
      V2 {},
      [](const V2& a, const V2& b) { return a + b; });
    st_.wmoved += sums[0];
    st_.facc += sums[1];

    if constexpr(std::same_as<Ch, ChSSA>) {
      const auto dgam = ptr(sam_.dgam_);

      // all sums are in; the first event to reach a lepton takes its sum and leaves
      // zero, so the lepton is rescaled once and dgam is clean for the next pass
      thrust::for_each(w.on_this(), ord0, ord0 + nev, [=](const I e) {
        const auto i  = rec[acc[e]].j;
        const auto dg = static_cast<double>(sstd::atomic_exchange(&dgam[i], 0.0f));
        if(dg <= 0.0 or not v.alive(i)) { return; }

        // rescale |u| to the energy g0 + dg, on the mass shell
        const auto u0 = v.u(i);
        const auto g0 = sstd::sqrt(1.0 + toolbox::dot(u0, u0));
        const auto g1 = g0 + dg;
        v.set_u(
          i,
          u0 *
            (sstd::sqrt(g1 * g1 - 1.0) / sstd::sqrt(sstd::max(g0 * g0 - 1.0, 1e-300))));

        v.touch(i);
      });
      w.wait();
    }

    if constexpr(std::same_as<Ch, ChCompton>) {
      const auto kick = ptr(sam_.ev_.kick);
      const auto n2p  = ptr(sam_.ev_.n2);

      // atomic: several events may kick one lepton
      thrust::for_each(w.on_this(), ord0, ord0 + nev, [=](const I e) {
        if(n2p[e] <= 0) { return; }
        const auto i = rec[acc[e]].j;
        for(auto k = 0uz; k < 3uz; ++k) {
          sstd::atomic_add(&v.vel[i][k], kick[3 * e + k]);
        }
      });

      // once all kicks are in; touching a lepton hit by several events repeats the same
      thrust::for_each(w.on_this(), ord0, ord0 + nev, [=](const I e) {
        const auto i = rec[acc[e]].j;
        if(v.alive(i)) { v.touch(i); }
      });
      w.wait();
    }
  }
};

void
  Sampler::twobody(Ctx& c, const Channel which)
{
  with_record(which, c.t->field.b_tile, box.x_min, [&](const auto& ch) {
    TwoBody { *this, c, ch, which }.run();
  });
}

// ---------------------------------------------------------------- diagnostics
double
  Sampler::pair_sum(const Channel which, const TileView& t, const bool weight) const
{
  const auto w  = tyvi::mdgrid_work {};
  const auto N2 = box.N2();
  const auto fv = t.field;  // the partner's field, for the absorption rate

  // brute-force sum over every (incident, partner) pair of the channel
  return with_record(which, fv.b_tile, box.x_min, [&](const auto& ch) {
    const auto i1 = static_cast<std::size_t>(std::remove_cvref_t<decltype(ch)>::t1);
    const auto i2 = static_cast<std::size_t>(ch.t2);
    auto& c1      = *t.con[i1];
    auto& c2      = *t.con[i2];
    const auto v1 = make_view(c1, box.w_lep);  // photons carry their own weight
    const auto v2 = make_view(c2, box.w_lep);
    const auto n1 = c1.size(), n2 = c2.size();
    const auto r0 = thrust::counting_iterator<I>(0);

    // one lane per incident slot, looping over every partner slot
    auto tot = thrust::transform_reduce(
      w.on_this(),
      r0,
      r0 + static_cast<std::ptrdiff_t>(n1),
      [=](const I i) {
        if(not v1.alive(i)) { return 0.0; }
        const auto p1 = v1.u(i);
        const auto wi = v1.w(i);
        auto s        = 0.0;

        for(I j = 0; j < n2; ++j) {
          if(not v2.alive(j)) { continue; }
          const auto wj = v2.w(j);
          s += sstd::max(wi, wj) * (weight ? sstd::min(wi, wj) : 1.0) * 2.0 *
               ch.rate(p1, v2.u(j), fv);
        }
        return s;
      },
      0.0,
      thrust::plus<double> {});

    if(i1 == i2) { tot *= 0.5; }  // each unordered pair once
    return tot / N2;
  });
}

}  // namespace qed
