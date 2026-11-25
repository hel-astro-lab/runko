// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

// The one-body channels: quantum synchrotron emission and magnetic pair production.

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

void
  Sampler::onebody(Ctx& c, const Channel which, const int l)
{
  const auto emit = which != Channel::bw;  // emission (syn) or pair production (bw)
  const auto is   = static_cast<std::size_t>(
    which == Channel::syn_m   ? Species::em
    : which == Channel::syn_p ? Species::ep
                              : Species::ph);

  auto& st     = stats.ch[static_cast<std::size_t>(l)];
  const auto w = tyvi::mdgrid_work {};

  c.refresh_all();
  const auto v = c.v[is];
  const auto n = c.t->con[is]->size();
  if(n == 0) { return; }

  const auto fv    = c.t->field;
  const auto N1    = box.N1();
  const auto x_min = box.x_min;
  if(not emit) {
    // no photon above the pair threshold: nothing to fire
    const auto lay  = v.lay;
    const auto k0   = thrust::counting_iterator<int>(0);
    const auto emax = thrust::transform_reduce(
      w.on_this(),
      k0,
      k0 + lay.nkeys,
      [=](const int k) { return lay.emax[k]; },
      0.0f,
      thrust::maximum<float> {});
    if(emax < static_cast<float>(bw::XMIN)) { return; }
  }

  // streams: A fires, C draws the products and their copies
  const auto keyA =
    runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(100 + 2 * l));
  const auto keyC =
    runko::rng_key(par.seed, istep_, 0, static_cast<std::uint64_t>(101 + 2 * l));

  st.M += static_cast<double>(n);

  const auto ds =
    thrust::raw_pointer_cast(dstat_.data()) + static_cast<std::size_t>(l) * d_count;
  staged_.resize(emit ? n : 2 * n);  // one photon per lepton; a pair per photon
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto cnt = thrust::raw_pointer_cast(counters_.data());
  thrust::fill(counters_.begin(), counters_.begin() + 1, 0);

  const auto a_res = par.a_res, x_ref = par.x_ref, acc_syn = par.acc_syn,
             kappa = par.kappa;
  const auto n_cap = par.n_cap;
  const auto wl    = box.w_lep;

  // append one product record at the parent's position and cell
  const auto stage = [=](
                       const Vec3& p,
                       const I i,
                       const double w_each,
                       const int n_r,
                       const int species,
                       const int bit) {
    const auto k = static_cast<std::size_t>(
      sstd::atomic_fetch_add(&cnt[c_staged], std::int32_t { 1 }));
    const float pos[3] = { v.pos[i][0], v.pos[i][1], v.pos[i][2] };
    stg[k]             = staged_at(
      p,
      pos,
      w_each,
      n_r,
      v.lay.cell_of(pos[0], pos[1], pos[2]),
      species,
      bit);
  };

  const auto s0 = thrust::counting_iterator<I>(0);
  if(emit) {
    // Resolved (below one emission per step): fire at the merged probability, emit one
    // photon of the merged weight. Unresolved: cool continuously, emit one photon LP
    // for all the emissions of the step. f_res then thins soft photons by copies().
    thrust::for_each(
      w.on_this(),
      s0,
      s0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return; }

        auto rngA      = runko::Rng { keyA, i };
        auto rngC      = runko::Rng { keyC, i };
        const auto u   = v.u(i);
        const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
        if(gam < syn::GMIN) { return; }  // below the channel's threshold

        const auto b     = fv.b_eff(v.pos[i][0], v.pos[i][1], v.pos[i][2], u / gam);
        const auto chi_e = gam * b;
        const auto s = syn::shape(chi_e);  // read by the rate, the mean and the sampler

        // emissions per step and mean photon energy, both above the floor
        const auto lam = syn::rate(s, chi_e, gam, x_min, false) / N1;
        if(lam <= 0.0) { return; }
        const auto xm = syn::xmean(s, gam, false, x_min);

        // slow: resolved; f: emissions merged into one LP, which fires with the
        // probability lam / f <= 1, so that the emissions of a step average lam
        const auto slow = lam <= 1.0;
        const auto f    = slow ? facc_syn(gam, xm, acc_syn, kappa) : 1.0;
        const auto fire = slow ? rngA.uniform() < lam / f : true;

        auto dgu = 0.0;  // energy the continuous cooling took
        if(not slow) {
          // the unresolved parent cools whether or not it emits an LP this step
          const auto g1 = sstd::max(gamma_cooled(gam, lam, xm, b, N1), 1.0);
          v.set_u(
            i,
            u * (sstd::sqrt(sstd::max(g1 * g1 - 1.0, 0.0)) /
                 sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
          dgu = gam - g1;
        }

        if(not fire) {
          if(not slow) { v.touch(i); }
          return;
        }

        // emissions the LP stands for; the photon is drawn from the pre-cooling state
        const auto mult = slow ? 1.0 : dgu / sstd::max(xm, 1e-300);
        const auto P    = syn::sample(s, chi_e, u, rngC, false, x_min);
        const auto x    = toolbox::norm(P.k);

        // the LP carries f x if resolved, the cooling loss if not; clipped at the
        // mass shell, the emitted weight shrinks with it
        const auto dg  = slow ? f * x : dgu;
        const auto dgc = sstd::min(dg, gam - 1.0);
        const auto wem =
          v.w(i) * f * mult * (dg > 0.0 ? dgc / sstd::max(dg, 1e-300) : 1.0);

        if(slow) {
          // a resolved parent recoils by what its photon carries
          const auto g1 = sstd::max(gam - dgc, 1.0);
          v.set_u(
            i,
            u * (sstd::sqrt(sstd::max(g1 * g1 - 1.0, 0.0)) /
                 sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
        }
        v.touch(i);

        sstd::atomic_add(&ds[d_fired], 1.0);
        sstd::atomic_add(&ds[d_wem], wem);
        sstd::atomic_add(&ds[d_facc], f);

        // the photon LP, split or rouletted by the resolution function
        const auto cp = copies(f_res(x, a_res, x_ref), wem, rngC, n_cap);
        if(cp.n_r >= 1) {
          stage(P.k, i, cp.w_each, cp.n_r, 0, static_cast<int>(rngC.uniform() < 0.5));
        }
      });
    w.wait();  // place() reads the staged count on the host
  } else {
    // magnetic pair production: the photon goes, pairs come at the lepton weight
    ev_q3_.resize(n);  // pairs to create per fired photon, w/w_lep
    ev_n3_.resize(n);  // its carry-rounded count
    ev_cell_.resize(n);
    lane_acc_.assign(n, 0);  // fired flag per slot
    const auto lacc = thrust::raw_pointer_cast(lane_acc_.data());

    // one lane per photon slot: fire on the conversion rate
    thrust::for_each(
      w.on_this(),
      s0,
      s0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return; }

        auto rngA    = runko::Rng { keyA, i };
        const auto k = v.u(i);
        const auto b = fv.b_eff(
          v.pos[i][0],
          v.pos[i][1],
          v.pos[i][2],
          k / sstd::max(toolbox::norm(k), 1e-300));

        const auto lam = bw::rate(k, b) / N1;  // conversions per step
        const auto z   = sstd::max(rngA.uniform(), 1e-300);
        if(lam > 0.0 and -sstd::log(z) < lam) { lacc[i] = 1; }
      });

    acc_.resize(n);
    const auto acc_end = thrust::copy_if(
      w.on_this(),
      s0,
      s0 + static_cast<std::ptrdiff_t>(n),
      acc_.begin(),
      [=](const I i) { return lacc[i] != 0; });

    const auto na = static_cast<std::size_t>(acc_end - acc_.begin());
    st.acc += static_cast<double>(na);
    if(na == 0) {
      w.wait();
      return;
    }

    const auto acc   = thrust::raw_pointer_cast(acc_.data());
    const auto q3    = thrust::raw_pointer_cast(ev_q3_.data());
    const auto n3    = thrust::raw_pointer_cast(ev_n3_.data());
    const auto ecell = thrust::raw_pointer_cast(ev_cell_.data());
    const auto e0    = thrust::counting_iterator<I>(0);

    // one lane per fired photon: the pair count it asks for and its cell
    thrust::for_each(
      w.on_this(),
      e0,
      e0 + static_cast<std::ptrdiff_t>(na),
      [=](const I e) {
        const auto i = acc[e];
        q3[e]        = v.w(i) / wl;
        ecell[e]     = v.lay.cell_of(v.pos[i][0], v.pos[i][1], v.pos[i][2]);
      });
    w.wait();

    // the pair counts, rounded systematically within every cell
    select(w, false, l, 2, c.ncell, na, ecell, q3, n3);

    // one lane per fired photon: draw the pair, stage its copies, remove the photon
    thrust::for_each(
      w.on_this(),
      e0,
      e0 + static_cast<std::ptrdiff_t>(na),
      [=](const I e) {
        const auto i = acc[e];
        auto rngC    = runko::Rng { keyC, i };
        const auto k = v.u(i);
        const auto b = fv.b_eff(
          v.pos[i][0],
          v.pos[i][1],
          v.pos[i][2],
          k / sstd::max(toolbox::norm(k), 1e-300));
        const auto P = bw::sample(k, b, rngC);

        if(n3[e] > 0) {
          // both charges at the photon's position: no charge density is added
          stage(P.um, i, wl, n3[e], 1, 0);
          stage(P.up, i, wl, n3[e], 2, 0);
        }
        v.kill(i);  // s = w: the parent always goes
      });
    w.wait();
  }

  this->place(c);
}

}  // namespace qed
