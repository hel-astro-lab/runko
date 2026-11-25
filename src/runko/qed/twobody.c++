// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

// The two-body channels: envelope, strata, thinning, reaction.

#include "runko/qed/sampler.h"

#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/qed/weights.h"
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

}  // namespace

template<typename Ch>
void
  Sampler::twobody(Ctx& c, const Ch& ch, const int l)
{
  // side 1 is the incident species, side 2 the partner
  const auto i1    = static_cast<std::size_t>(Ch::t1);
  const auto i2    = static_cast<std::size_t>(ch.t2);
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;  // (cell, bin) groups
  const auto npair = ncb * nb;    // (cell, bin a, bin b) pairs

  auto& st = stats.ch[static_cast<std::size_t>(l)];
  const auto ds =
    thrust::raw_pointer_cast(dstat_.data()) + static_cast<std::size_t>(l) * d_count;

  const auto w = tyvi::mdgrid_work {};
  c.refresh_all();

  // ---- per side and (cell, bin): the slot sub-ranges and the bounds over them
  // base0, cnt0: [0, ncb) the sub-range of the first key, [ncb, 2 ncb) of the second
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
    const auto b0  = thrust::raw_pointer_cast(base0.data()),
               c0  = thrust::raw_pointer_cast(cnt0.data());
    const auto wm  = thrust::raw_pointer_cast(wmax.data());
    const auto lo  = thrust::raw_pointer_cast(emin.data()),
               hi  = thrust::raw_pointer_cast(emax.data());
    const auto inf = std::numeric_limits<float>::infinity();

    thrust::for_each(
      w.on_this(),
      thrust::counting_iterator<int>(0),
      thrust::counting_iterator<int>(ncb),
      [=](const int cb) {
        int k0 = -1, k1 = -1;  // the group's one or two keys
        switch(sel) {
          case Sel::cell: k0 = cb; break;
          case Sel::bit0: k0 = 2 * cb; break;
          case Sel::bit1: k0 = 2 * cb + 1; break;
          case Sel::both:
            k0 = 2 * cb;
            k1 = 2 * cb + 1;
            break;
        }

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
  };

  side(i1, Ch::sel1, base1_, cnt1_, wmax1_, emin1_, emax1_);
  side(i2, Ch::sel2, base2_, cnt2_, wmax2_, emin2_, emax2_);

  const auto s1 = SideView { thrust::raw_pointer_cast(base1_.data()),
                             thrust::raw_pointer_cast(cnt1_.data()),
                             thrust::raw_pointer_cast(base1_.data()) + ncb,
                             thrust::raw_pointer_cast(cnt1_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax1_.data()),
                             thrust::raw_pointer_cast(emin1_.data()),
                             thrust::raw_pointer_cast(emax1_.data()) };

  const auto s2 = SideView { thrust::raw_pointer_cast(base2_.data()),
                             thrust::raw_pointer_cast(cnt2_.data()),
                             thrust::raw_pointer_cast(base2_.data()) + ncb,
                             thrust::raw_pointer_cast(cnt2_.data()) + ncb,
                             thrust::raw_pointer_cast(wmax2_.data()),
                             thrust::raw_pointer_cast(emin2_.data()),
                             thrust::raw_pointer_cast(emax2_.data()) };

  // ---- envelope per bin pair, from the slot counts and the bins' bounds alone
  // per bin pair: envelope Lam, rate bound sh, weight bound wmp; Lam summed over b
  // into rowsum, over a into colsum
  Lam_.resize(static_cast<std::size_t>(npair));
  sh_.resize(static_cast<std::size_t>(npair));
  wm_.resize(static_cast<std::size_t>(npair));
  rowsum_.resize(static_cast<std::size_t>(ncb));
  colsum_.resize(static_cast<std::size_t>(ncb));

  const auto Lam    = thrust::raw_pointer_cast(Lam_.data()),
             sh     = thrust::raw_pointer_cast(sh_.data());
  const auto wmp    = thrust::raw_pointer_cast(wm_.data());
  const auto rowsum = thrust::raw_pointer_cast(rowsum_.data()),
             colsum = thrust::raw_pointer_cast(colsum_.data());

  const auto kappa = par.kappa;
  const auto N2c   = box.N2() / ncell;  // the two-body clock of one reaction cell

  thrust::for_each(
    w.on_this(),
    thrust::counting_iterator<int>(0),
    thrust::counting_iterator<int>(ncb),
    [=](const int ca) {
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
        if constexpr(Ch::has_dgm) {
          const auto glo = sstd::max(lo2, Ch::g_min);
          fk             = sstd::min(
            fkin(kappa, glo, ch.dgmax(hi1, glo)),
            fkin(kappa, hi2, ch.dgmax(hi1, hi2)));
        }
        const auto wmv = pair_weight(w1, w2, fk);

        // envelope: all n1 n2 slot pairs at the bin pair's rate and weight bounds
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
      sstd::atomic_add(&ds[d_shnum], sum_sh);
      sstd::atomic_add(&ds[d_shden], sum);
    });

  thrust::for_each(
    w.on_this(),
    thrust::counting_iterator<int>(0),
    thrust::counting_iterator<int>(ncb),
    [=](const int cb) {
      const auto cc = cb / nb, bb = cb % nb;
      auto sum = 0.0;
      for(auto a = 0; a < nb; ++a) { sum += Lam[(cc * nb + a) * nb + bb]; }
      colsum[cb] = sum;
    });

  stats.cell_lanes += npair;
  const auto g0 = thrust::counting_iterator<int>(0);

  // ---- sub-steps: distinct strata need at most as many events per bin as slots,
  // on both sides; beyond that the channel splits into sub-steps of a larger N_2.
  // Additive channels cap them: past the cap the row clamp drops the surplus (an
  // under-count the cap controls), so no two lanes ever hold one slot.
  auto nsub = thrust::transform_reduce(
    w.on_this(),
    g0,
    g0 + ncb,
    [=](const int ca) {
      const auto n = static_cast<double>(s1.N(ca));
      return rowsum[ca] > 0.0 and n >= 1.0
               ? sstd::ceil(rowsum[ca] / sstd::max(n - 1.0, 1.0))
               : 1.0;
    },
    1.0,
    thrust::maximum<double> {});

  if constexpr(Ch::kind == Kind::convert) {
    // a consumed partner needs distinct strata too
    nsub = sstd::max(
      nsub,
      thrust::transform_reduce(
        w.on_this(),
        g0,
        g0 + ncb,
        [=](const int cb) {
          const auto n = static_cast<double>(s2.N(cb));
          return colsum[cb] > 0.0 and n >= 1.0
                   ? sstd::ceil(colsum[cb] / sstd::max(n - 1.0, 1.0))
                   : 1.0;
        },
        1.0,
        thrust::maximum<double> {}));
  } else {
    nsub = sstd::min(nsub, static_cast<double>(par.nsub_max));
  }

  const auto ns    = static_cast<int>(nsub);
  st.nsub_max_seen = sstd::max(st.nsub_max_seen, ns);
  st.nsub_steps += (ns > 1);

  // each sub-step with the envelope divided by ns
  for(auto r = 0; r < ns; ++r) { pass(c, ch, l, r, ns, s1, s2); }
}

template<typename Ch>
void
  Sampler::pass(
    Ctx& c,
    const Ch& ch,
    const int l,
    const int r,
    const int nsub,
    const SideView& s1,
    const SideView& s2)
{
  const auto i1    = static_cast<std::size_t>(Ch::t1);
  const auto i2    = static_cast<std::size_t>(ch.t2);
  const auto nb    = c.nb;
  const auto ncell = c.ncell;
  const auto ncb   = ncell * nb;
  const auto npair = ncb * nb;

  auto& st        = stats.ch[static_cast<std::size_t>(l)];
  const auto w    = tyvi::mdgrid_work {};
  const auto step = istep_;

  // random streams of this (channel, sub-step): A the counts (lane = incident group)
  // and the shifts (lanes past them), B the thinning, C the reactions
  const auto keyA = runko::rng_key(
    par.seed,
    step,
    static_cast<std::uint64_t>(r),
    static_cast<std::uint64_t>(1 + 4 * l));
  const auto keyB = runko::rng_key(
    par.seed,
    step,
    static_cast<std::uint64_t>(r),
    static_cast<std::uint64_t>(2 + 4 * l));
  const auto keyC = runko::rng_key(
    par.seed,
    step,
    static_cast<std::uint64_t>(r),
    static_cast<std::uint64_t>(3 + 4 * l));

  const auto ds =
    thrust::raw_pointer_cast(dstat_.data()) + static_cast<std::size_t>(l) * d_count;

  const auto Lam = thrust::raw_pointer_cast(Lam_.data()),
             sh  = thrust::raw_pointer_cast(sh_.data());
  const auto wmp = thrust::raw_pointer_cast(wm_.data());

  // ---- counts per bin pair, clamps, strata bases and shifts
  // M events per bin pair; for a consumed partner off2, its strata base along a, and
  // M2, the events per partner bin
  M_.resize(static_cast<std::size_t>(npair));
  off2_.resize(static_cast<std::size_t>(npair));
  M2_.resize(static_cast<std::size_t>(ncb));

  const auto M    = thrust::raw_pointer_cast(M_.data());
  const auto off2 = thrust::raw_pointer_cast(off2_.data());
  const auto M2   = thrust::raw_pointer_cast(M2_.data());
  const auto cb0  = thrust::counting_iterator<int>(0);

  // incident rows: M[a, b] by cumulative rounding with one uniform for the whole row
  thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int ca) {
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

    {
      // row clamp: no more events than slots, so the strata stay distinct; a
      // backstop the sub-steps make rare, counted as clamped
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

  // partner columns of a consumed partner: the clamp, then the strata bases along a
  if constexpr(Ch::kind == Kind::convert) {
    thrust::for_each(w.on_this(), cb0, cb0 + ncb, [=](const int cb) {
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
        sstd::atomic_add(&ds[d_clamped], 1.0);
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

  // first lane of every bin pair; the last entry, the lane total, is read on the host.
  // The lanes run row by row, so a lane's rank in its row is its incident stratum.
  Mstart_.resize(static_cast<std::size_t>(npair) + 1);
  Mstart_[0] = I { 0 };
  thrust::inclusive_scan(w.on_this(), M_.begin(), M_.end(), Mstart_.begin() + 1);
  w.wait();

  const auto Mtot = static_cast<std::size_t>(Mstart_[static_cast<std::size_t>(npair)]);
  st.M += static_cast<double>(Mtot);

  if(Mtot == 0) { return; }

  // ---- thinning kernel: one lane per envelope event
  // per lane: the bin pair, the state, the incident and partner slot, the
  // accumulation factor
  lane_acc_.resize(Mtot);
  lane_i_.resize(Mtot);
  lane_j_.resize(Mtot);
  lane_S_.resize(Mtot);
  lane_g_.resize(Mtot);

  const auto Mstart = thrust::raw_pointer_cast(Mstart_.data());
  const auto lg     = thrust::raw_pointer_cast(lane_g_.data());
  const auto lacc   = thrust::raw_pointer_cast(lane_acc_.data());
  const auto li     = thrust::raw_pointer_cast(lane_i_.data()),
             lj     = thrust::raw_pointer_cast(lane_j_.data());
  const auto lS     = thrust::raw_pointer_cast(lane_S_.data());

  const auto v1 = c.v[i1], v2 = c.v[i2];
  const auto fv    = c.t->field;  // the field at the partner, for the absorption rate
  const auto kappa = par.kappa;
  const auto lane0 = thrust::counting_iterator<I>(0);

  thrust::for_each(
    w.on_this(),
    lane0,
    lane0 + static_cast<std::ptrdiff_t>(Mtot),
    [=](const I lane) {
      auto rng      = runko::Rng { keyB, lane };
      const auto g  = segment_of(Mstart, npair, lane);
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
      const auto row = Mstart[cb1 * nb];
      const auto p1  = static_cast<std::uint64_t>(lane - row);
      const auto M1  = static_cast<std::uint64_t>(Mstart[(cb1 + 1) * nb] - row);
      const auto vv1 = (shift(1, cb1, n1) + stratum(p1, M1, n1, rng)) % n1;

      I vv2;
      if constexpr(Ch::kind == Kind::convert) {
        // a consumed partner: a slot drawn in stratum p2 of the column's M2
        const auto p2 = static_cast<std::uint64_t>(off2[g]) +
                        static_cast<std::uint64_t>(lane - Mstart[g]);
        const auto M2c = static_cast<std::uint64_t>(M2[cb2]);
        vv2            = (shift(0, cb2, n2) + stratum(p2, M2c, n2, rng)) % n2;
      } else {
        // the partner only accumulates what the events give it, so several lanes may
        // share it: drawn uniformly with replacement
        vv2 = static_cast<I>(rng.uniform() * static_cast<double>(n2)) % n2;
      }

      const auto i = s1.pick(cb1, vv1), j = s2.pick(cb2, vv2);
      li[lane]   = i;
      lj[lane]   = j;
      lg[lane]   = g;
      lS[lane]   = 0.0;
      lacc[lane] = L_WFAIL;

      // a dead slot, or one a kick dropped from its key: nothing to react
      if(not(v1.alive(i) and v2.alive(j)) or v1.lay.skey[i] < 0 or v2.lay.skey[j] < 0) {
        lacc[lane] = L_VIRT;
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
      lacc[lane] = L_WPASS;

      const auto p1v = v1.u(i), p2v = v2.u(j);
      const auto rate = pair_rate(ch, fv, p1v, p2v);
      auto fk         = 1.0;
      if constexpr(Ch::has_dgm) {
        const auto x   = toolbox::norm(p1v);
        const auto gam = sstd::sqrt(1.0 + toolbox::dot(p2v, p2v));
        fk             = fkin(kappa, gam, ch.dgmax(x, gam));
      }

      lS[lane] = event_size(wi, wj, fk);
      if(u < (pair_weight(wi, wj, fk) / wm) * rate / sh[g]) { lacc[lane] = L_ACC; }
    });

  acc_.resize(Mtot);
  const auto acc_end = thrust::copy_if(
    w.on_this(),
    lane0,
    lane0 + static_cast<std::ptrdiff_t>(Mtot),
    acc_.begin(),
    [=](const I lane) { return lacc[lane] == L_ACC; });
  const auto na = static_cast<std::size_t>(acc_end - acc_.begin());  // accepted events
  st.acc += static_cast<double>(na);

  // the dead lanes and those past the weight stage only, in one reduction
  using V2         = toolbox::VecD<double, 2>;
  const auto vsum  = [](const V2& a, const V2& b) { return a + b; };
  const auto lanes = thrust::transform_reduce(
    w.on_this(),
    lane0,
    lane0 + static_cast<std::ptrdiff_t>(Mtot),
    [=](const I lane) {
      return V2(lacc[lane] == L_VIRT ? 1.0 : 0.0, lacc[lane] == L_WPASS ? 1.0 : 0.0);
    },
    V2 {},
    vsum);

  st.virt += lanes[0];
  st.wpass += static_cast<double>(na) + lanes[1];

  if(na == 0) {
    w.wait();
    return;
  }

  // ---- decide: what an event does to its parents follows from their weights and its
  // size alone, so it is settled before the final state is drawn
  // q1, q2: the update or removal probability of parent 1, 2; q3: the lepton pairs
  // to create (photon-photon); n1p, n2p, n3p: their carry-rounded decisions; ef: the
  // exchanges the partner takes as one
  ev_q1_.resize(na);
  ev_q2_.resize(na);
  ev_q3_.resize(na);
  ev_n1_.resize(na);
  ev_n2_.resize(na);
  ev_n3_.resize(na);
  ev_cell_.resize(na);
  ev_f_.resize(na);

  const auto acc   = thrust::raw_pointer_cast(acc_.data());
  const auto q1    = thrust::raw_pointer_cast(ev_q1_.data()),
             q2    = thrust::raw_pointer_cast(ev_q2_.data());
  const auto q3    = thrust::raw_pointer_cast(ev_q3_.data());
  const auto n1p   = thrust::raw_pointer_cast(ev_n1_.data()),
             n2p   = thrust::raw_pointer_cast(ev_n2_.data());
  const auto n3p   = thrust::raw_pointer_cast(ev_n3_.data());
  const auto ecell = thrust::raw_pointer_cast(ev_cell_.data());
  const auto ef    = thrust::raw_pointer_cast(ev_f_.data());

  const auto wl  = box.w_lep;
  const auto ev0 = thrust::counting_iterator<I>(0);

  thrust::for_each(
    w.on_this(),
    ev0,
    ev0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) {
      const auto lane = acc[e];
      const auto wi = v1.w(li[lane]), wj = v2.w(lj[lane]);
      const auto S  = lS[lane];  // the weight the event moves, <= wi
      const auto n2 = S / wj;    // the exchanges the partner takes
      ecell[e]      = lg[lane] / (nb * nb);

      // the incident gives up S; the partner takes n2 exchanges, with that probability
      // while n2 < 1 and as one event beyond
      q1[e] = S / wi;
      q2[e] = sstd::min(n2, 1.0);
      ef[e] = sstd::max(n2, 1.0);
      q3[e] = Ch::pairs ? S / wl : 0.0;
    });
  // no wait: the rounding and the reactions follow on the same queue

  // the systematic rounding of the probabilities in every cell
  select(w, true, l, 0, ncell, na, ecell, q1, n1p);
  select(w, true, l, 1, ncell, na, ecell, q2, n2p);
  if constexpr(Ch::pairs) { select(w, true, l, 2, ncell, na, ecell, q3, n3p); }

  // ---- react: the final state of every event; the incident is updated in place, the
  // products are staged. Lanes may share a partner, so its state is only read here:
  // what the events give it is summed (kicks per event, absorbed energy per slot) and
  // applied once the kernel is through.
  staged_.resize(2 * na);  // at most two products per event
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());

  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);

  if constexpr(Ch::kind == Kind::scatter) { ev_kick_.resize(3 * na); }

  // dgam is zero between passes (apply_dgam clears what it reads): only growth needs
  // new zeros
  if constexpr(Ch::kind == Kind::absorb) {
    if(dgam_.size() < c.t->con[i2]->size()) {
      dgam_.resize(c.t->con[i2]->size(), 0.0f);
    }
  }

  const auto kick  = thrust::raw_pointer_cast(ev_kick_.data());
  const auto dgam  = thrust::raw_pointer_cast(dgam_.data());
  const auto a_res = par.a_res, x_ref = par.x_ref;
  const auto n_cap = par.n_cap;

  // append one product record to the staging buffer
  const auto stage = [=](
                       const Vec3& p,
                       const float* pos,
                       const double w_each,
                       const int n_r,
                       const int cell,
                       const int species,
                       const int bit) {
    const auto k = static_cast<std::size_t>(
      sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    stg[k] = staged_at(p, pos, w_each, n_r, cell, species, bit);
  };

  thrust::for_each(
    w.on_this(),
    ev0,
    ev0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) {
      auto rng        = runko::Rng { keyC, e };
      const auto lane = acc[e];
      const auto i = li[lane], j = lj[lane];
      const auto f    = ef[e];
      const auto cell = ecell[e];
      const auto p1v = v1.u(i), p2v = v2.u(j);

      // products are born where the incident is
      const float posi[3] = { v1.pos[i][0], v1.pos[i][1], v1.pos[i][2] };

      [[maybe_unused]] auto ntry = 0;  // proposals of the rejection samplers

      if constexpr(Ch::kind == Kind::scatter) {
        // Compton: the lepton takes the f exchanges as one momentum kick
        const auto P  = comp::sample(p1v, p2v, rng, ntry);
        const auto dp = (-f) * (P.k1 - p1v);
        for(auto k = 0uz; k < 3uz; ++k) { kick[3 * e + k] = static_cast<float>(dp[k]); }

        if(n1p[e] > 0) {
          // the parent takes the new state as the copies the resolution function asks
          // for, itself the first copy
          const auto x0 = v1.e(i), x1 = toolbox::norm(P.k1);
          const auto n  = f_res(x1, a_res, x_ref) / f_res(x0, a_res, x_ref);
          const auto cp = copies(n, v1.w(i), rng, n_cap);

          if(cp.n_r >= 1) {
            // read before the touch, which drops a parent that left its bin from its key
            const auto bit = v1.lay.skey[i] % 2;
            v1.set_u(i, P.k1);
            v1.set_w(i, cp.w_each);
            v1.touch(i);
            if(cp.n_r > 1) { stage(P.k1, posi, cp.w_each, cp.n_r - 1, cell, 0, bit); }
          } else {
            v1.kill(i);  // the roulette merged it away
          }
        }
      } else if constexpr(Ch::kind == Kind::absorb) {
        // self-absorption: the lepton gains the energy f x, summed per lepton
        if(n2p[e] > 0) {
          sstd::atomic_add(&dgam[j], static_cast<float>(f * toolbox::norm(p1v)));
        }
        if(n1p[e] > 0) { v1.kill(i); }
      } else if constexpr(Ch::pairs) {
        // photon-photon: the weight s becomes lepton pairs at the lepton weight
        const auto P    = gg::sample(p1v, p2v, rng, ntry);
        const auto flip = rng.uniform() < 0.5;
        const auto nl   = n3p[e];
        if(nl > 0) {
          // a pair born on top of itself adds no charge density: div E is unchanged
          stage(flip ? P.p1 : P.p2, posi, wl, nl, cell, 1, 0);
          stage(flip ? P.p2 : P.p1, posi, wl, nl, cell, 2, 0);
        }
        if(n1p[e] > 0) { v1.kill(i); }
        if(n2p[e] > 0) { v2.kill(j); }
      } else {
        // pair annihilation: two photons of weight s, each realized as the copies its
        // resolution asks for
        const auto P = ann::sample(p1v, p2v, rng, ntry);
        const auto s = sstd::min(v1.w(i), v2.w(j));
        auto photon  = [&](const Vec3& k) {
          const auto cp = copies(f_res(toolbox::norm(k), a_res, x_ref), s, rng, n_cap);
          const auto bit = static_cast<int>(rng.uniform() < 0.5);
          if(cp.n_r >= 1) { stage(k, posi, cp.w_each, cp.n_r, cell, 0, bit); }
        };
        photon(P.k0);
        photon(P.k1);
        if(n1p[e] > 0) { v1.kill(i); }
        if(n2p[e] > 0) { v2.kill(j); }
      }
    });
  w.wait();

  const auto sums = thrust::transform_reduce(
    w.on_this(),
    ev0,
    ev0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) { return V2(lS[acc[e]], ef[e]); },
    V2 {},
    vsum);
  st.wmoved += sums[0];
  st.facc += sums[1];

  // ---- the partners take their summed changes, then the products are placed
  if constexpr(Ch::kind == Kind::absorb) { apply_dgam(c, ch.t2, na); }
  if constexpr(Ch::kind == Kind::scatter) { apply_kicks(c, ch.t2, na); }
  this->place(c);
}

void
  Sampler::twobody(Ctx& c, const Channel which)
{
  const auto l      = static_cast<int>(which);
  const auto b_tile = c.t->field.b_tile;

  switch(which) {
    case Channel::compton_m: twobody(c, ChCompton { Species::em }, l); break;
    case Channel::compton_p: twobody(c, ChCompton { Species::ep }, l); break;
    case Channel::pairann: twobody(c, ChAnn {}, l); break;
    case Channel::photann: twobody(c, ChGG {}, l); break;
    case Channel::ssa_m:
      twobody(c, ChSSA { Species::em, b_tile, box.x_min }, l);
      break;
    case Channel::ssa_p:
      twobody(c, ChSSA { Species::ep, b_tile, box.x_min }, l);
      break;
    default: throw std::runtime_error { "twobody: not a two-body channel" };
  }
}

void
  Sampler::apply_dgam(Ctx& c, const Species s, const std::size_t na)
{
  const auto v    = c.v[static_cast<std::size_t>(s)];
  const auto dgam = thrust::raw_pointer_cast(dgam_.data());
  const auto acc  = thrust::raw_pointer_cast(acc_.data());
  const auto lj   = thrust::raw_pointer_cast(lane_j_.data());
  const auto w    = tyvi::mdgrid_work {};
  const auto e0   = thrust::counting_iterator<I>(0);

  // all sums are in; the first event to reach a lepton takes its sum and leaves zero,
  // so the lepton is rescaled once and dgam is clean for the next pass
  thrust::for_each(
    w.on_this(),
    e0,
    e0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) {
      const auto i  = lj[acc[e]];
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

void
  Sampler::apply_kicks(Ctx& c, const Species s, const std::size_t na)
{
  const auto v    = c.v[static_cast<std::size_t>(s)];
  const auto acc  = thrust::raw_pointer_cast(acc_.data());
  const auto lj   = thrust::raw_pointer_cast(lane_j_.data());
  const auto kick = thrust::raw_pointer_cast(ev_kick_.data());
  const auto n2p  = thrust::raw_pointer_cast(ev_n2_.data());
  const auto w    = tyvi::mdgrid_work {};
  const auto e0   = thrust::counting_iterator<I>(0);

  // atomic: several events may kick one lepton
  thrust::for_each(
    w.on_this(),
    e0,
    e0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) {
      if(n2p[e] <= 0) { return; }
      const auto i = lj[acc[e]];
      for(auto k = 0uz; k < 3uz; ++k) { sstd::atomic_add(&v.vel[i][k], kick[3 * e + k]); }
    });

  // once all kicks are in; touching a lepton hit by several events repeats the same
  thrust::for_each(
    w.on_this(),
    e0,
    e0 + static_cast<std::ptrdiff_t>(na),
    [=](const I e) {
      const auto i = lj[acc[e]];
      if(v.alive(i)) { v.touch(i); }
    });
  w.wait();
}


// ---------------------------------------------------------------- diagnostics
double
  Sampler::pair_sum(const Channel which, const TileView& t, const bool weight) const
{
  const auto w  = tyvi::mdgrid_work {};
  const auto N2 = box.N2();
  const auto fv = t.field;  // the partner's field, for the absorption rate

  // brute-force sum over every (incident, partner) pair of the channel's record type
  auto run = [&](const auto& ch) {
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
               pair_rate(ch, fv, p1, v2.u(j));
        }
        return s;
      },
      0.0,
      thrust::plus<double> {});

    if(i1 == i2) { tot *= 0.5; }  // each unordered pair once
    return tot / N2;
  };

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

}  // namespace qed
