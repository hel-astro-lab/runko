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

namespace {

/// Stages a product at the position and in the reaction cell of its parent.
struct [[nodiscard]] AtParent {
  PView v;  // the parent species
  Stager out;

  constexpr void operator()(
    const Vec3& p,
    const I i,
    const double w_each,
    const int n_r,
    const Species s,
    const int bit) const
  {
    const float pos[3] = { v.pos[i][0], v.pos[i][1], v.pos[i][2] };
    out(p, pos, w_each, n_r, v.lay.cell_of(pos[0], pos[1], pos[2]), s, bit);
  }
};

}  // namespace

// ---------------------------------------------------------------- synchrotron emission
void
  Sampler::synchrotron(Ctx& c, const Channel which)
{
  const auto l  = static_cast<std::size_t>(which);
  const auto is = static_cast<std::size_t>(
    which == Channel::syn_m ? Species::em : Species::ep);

  c.refresh_all();
  const auto v = c.v[is];
  const auto n = c.t->con[is]->size();
  if(n == 0) { return; }

  stats.ch[l].M += static_cast<double>(n);

  // random streams, the lane being the slot: A fires, C draws the photon and its copies
  const auto keyA  = runko::rng_key(par.seed, istep_, 0, stream::onebody(which, 0));
  const auto keyC  = runko::rng_key(par.seed, istep_, 0, stream::onebody(which, 1));
  const auto ds    = ptr(dstat_) + l;
  const auto stage = AtParent { v, this->stager(n) };  // one photon per lepton

  const auto fv    = c.t->field;
  const auto N1    = box.N1();
  const auto x_min = box.x_min;
  const auto a_res = par.a_res, x_ref = par.x_ref, acc_syn = par.acc_syn,
             kappa = par.kappa;
  const auto n_cap = par.n_cap;

  const auto w  = tyvi::mdgrid_work {};
  const auto s0 = thrust::counting_iterator<I>(0);

  // one lane per lepton slot: plan the step, draw the photon, recoil, stage the LP
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
      if(gam < syn::gamma_min) { return; }  // below the channel's threshold

      const auto b     = fv.b_eff(v.pos[i][0], v.pos[i][1], v.pos[i][2], u / gam);
      const auto chi_e = gam * b;
      const auto s = syn::shape(chi_e);  // read by the rate, the mean and the sampler

      // emissions per step and mean photon energy, both above the floor
      const auto lam = syn::rate(s, chi_e, gam, x_min, false) / N1;
      if(lam <= 0.0) { return; }
      const auto xm = syn::xmean(s, gam, false, x_min);

      const auto plan =
        emission_plan(gam, lam, xm, b, N1, acc_syn, kappa, rngA.uniform());
      if(not plan.fire) { return; }

      // the photon is drawn from the state the step began with
      const auto P  = syn::sample(s, chi_e, u, rngC, false, x_min);
      const auto x  = toolbox::norm(P.k);
      const auto em = emitted(plan, v.w(i), gam, xm, x);

      // the parent keeps its direction
      v.set_u(
        i,
        u * (sstd::sqrt(sstd::max(em.g1 * em.g1 - 1.0, 0.0)) /
             sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300))));
      v.touch(i);

      sstd::atomic_add(&ds->fired, 1.0);
      sstd::atomic_add(&ds->wem, em.wem);
      sstd::atomic_add(&ds->facc, plan.f);

      // the photon LP, split or rouletted by the resolution function
      const auto cp = copies(f_res(x, a_res, x_ref), em.wem, rngC, n_cap);
      if(cp.n_r >= 1) {
        const auto bit = static_cast<int>(rngC.uniform() < 0.5);
        stage(P.k, i, cp.w_each, cp.n_r, Species::ph, bit);
      }
    });
  w.wait();  // place() reads the staged count on the host

  this->place(c);
}

// ---------------------------------------------------------------- magnetic pair production
void
  Sampler::breit_wheeler(Ctx& c)
{
  const auto l  = static_cast<std::size_t>(Channel::bw);
  const auto is = static_cast<std::size_t>(Species::ph);

  c.refresh_all();
  const auto v = c.v[is];
  const auto n = c.t->con[is]->size();
  if(n == 0) { return; }

  const auto w = tyvi::mdgrid_work {};
  {
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
    if(emax < static_cast<float>(bw::x_threshold)) { return; }
  }

  auto& st = stats.ch[l];
  st.M += static_cast<double>(n);

  // random streams, the lane being the slot: A fires, C draws the pair
  const auto keyA  = runko::rng_key(par.seed, istep_, 0, stream::onebody(Channel::bw, 0));
  const auto keyC  = runko::rng_key(par.seed, istep_, 0, stream::onebody(Channel::bw, 1));
  const auto stage = AtParent { v, this->stager(2 * n) };  // a pair per photon

  const auto fv = c.t->field;
  const auto N1 = box.N1();
  const auto wl = box.w_lep;

  ev_.q3.resize(n);  // pairs to create per fired photon, w/w_lep
  ev_.n3.resize(n);  // its carry-rounded count
  ev_.cell.resize(n);
  fired_.assign(n, 0);
  const auto fired = ptr(fired_);
  const auto s0    = thrust::counting_iterator<I>(0);
  const auto nslot = static_cast<std::ptrdiff_t>(n);

  // one lane per photon slot: the photon converts within the step with the
  // probability 1 - exp(-lam), exact for a decay
  thrust::for_each(w.on_this(), s0, s0 + nslot, [=](const I i) {
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
    if(lam > 0.0 and -sstd::log(z) < lam) { fired[i] = 1; }
  });

  acc_.resize(n);
  const auto acc_end =
    thrust::copy_if(w.on_this(), s0, s0 + nslot, acc_.begin(), [=](const I i) {
      return fired[i] != 0;
    });

  const auto na = static_cast<std::size_t>(acc_end - acc_.begin());
  st.acc += static_cast<double>(na);
  if(na == 0) {
    w.wait();
    return;
  }

  const auto acc   = ptr(acc_);
  const auto q3    = ptr(ev_.q3);
  const auto n3    = ptr(ev_.n3);
  const auto ecell = ptr(ev_.cell);
  const auto nev   = static_cast<std::ptrdiff_t>(na);

  // one lane per fired photon: the pair count it asks for and its cell
  thrust::for_each(w.on_this(), s0, s0 + nev, [=](const I e) {
    const auto i = acc[e];
    q3[e]        = v.w(i) / wl;
    ecell[e]     = v.lay.cell_of(v.pos[i][0], v.pos[i][1], v.pos[i][2]);
  });
  w.wait();

  // the pair counts, rounded systematically within every cell
  select(w, false, static_cast<int>(l), 2, c.ncell, na, ecell, q3, n3);

  // one lane per fired photon: draw the pair, stage its copies, remove the photon
  thrust::for_each(w.on_this(), s0, s0 + nev, [=](const I e) {
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
      stage(P.um, i, wl, n3[e], Species::em, 0);
      stage(P.up, i, wl, n3[e], Species::ep, 0);
    }
    v.kill(i);  // the whole weight converts: the parent always goes
  });
  w.wait();

  this->place(c);
}

}  // namespace qed
