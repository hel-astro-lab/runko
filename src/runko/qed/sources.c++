// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

// Sources and sinks of the open box: injection, heating, escape and their budget.

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
  Sampler::inject_pairs_plaw(
    TileView& t,
    const double ell,
    const double p,
    const double gmin,
    const double gmax)
{
  this->lay_if_needed(t);

  // pairs per step = energy the compactness injects / (2 <gamma> w_lep), rounded
  // stochastically
  double gm;  // mean gamma of the power law
  if(p == 1.0) {
    gm = (gmax - gmin) / std::log(gmax / gmin);
  } else {
    const auto a = 1.0 - p, b = 2.0 - p;
    gm = p == 2.0 ? (b / a) * (std::pow(gmax, a) - std::pow(gmin, a)) /
                      (std::pow(gmax, b) - std::pow(gmin, b))
                  : (a / b) * (std::pow(gmax, b) - std::pow(gmin, b)) /
                      (std::pow(gmax, a) - std::pow(gmin, a));
  }

  const auto wl    = box.w_lep;
  const auto npair = box.dE(ell) / (2.0 * gm * wl);
  auto rng0        = runko::Rng { runko::rng_key(par.seed, istep_, 0, stream::pairs_count), 0 };
  const auto n     = static_cast<std::size_t>(std::floor(npair + rng0.uniform()));
  if(n == 0) { return; }

  Ctx c { t, par.nbin };
  const auto stage = this->stager(2 * n);
  tally_[0] = Tally {};

  const auto db   = ptr(tally_);
  const auto key  = runko::rng_key(par.seed, istep_, 0, stream::pairs);
  const auto mins = t.mins, maxs = t.maxs;
  const auto lay = c.v[1].lay;
  const auto w   = tyvi::mdgrid_work {};
  const auto i0  = thrust::counting_iterator<I>(0);

  // the cumulative of gamma^-p at the ends
  const auto glo = sstd::pow(gmin, 1.0 - p), ghi = sstd::pow(gmax, 1.0 - p);

  // one lane per pair: uniform position; the two charges draw power-law energies and
  // isotropic directions independently
  thrust::for_each(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    [=](const I i) {
      auto rng           = runko::Rng { key, i };
      const float pos[3] = {
        mins[0] + (maxs[0] - mins[0]) * static_cast<float>(rng.uniform()),
        mins[1] + (maxs[1] - mins[1]) * static_cast<float>(rng.uniform()),
        mins[2] + (maxs[2] - mins[2]) * static_cast<float>(rng.uniform())
      };

      const auto cell = lay.cell_of(pos[0], pos[1], pos[2]);

      // electron and positron at the same position: neutral
      for(auto sp = 1; sp <= 2; ++sp) {
        const auto z   = rng.uniform();
        const auto gam = sstd::pow(
          glo + z * (ghi - glo),
          1.0 / (1.0 - p));  // inverse cumulative of gamma^-p
        const auto u = sstd::sqrt(gam * gam - 1.0) * iso_dir(rng);
        stage(u, pos, wl, 1, cell, static_cast<Species>(sp), 0);
        sstd::atomic_add(&db->E, wl * gam);
      }
    });
  w.wait();

  {
    const Tally h = tally_[0];
    budget.inj_E += h.E;
    budget.inj_W += 2.0 * static_cast<double>(n) * wl;
  }

  this->place(c);
}

void
  Sampler::inject_photons(
    TileView& t,
    const double ell,
    const double theta,
    const double wph,
    const bool planck)
{
  this->lay_if_needed(t);

  // photons per step = energy the compactness injects / (<x> wph), rounded
  // stochastically
  const auto xmean = (planck ? 2.701 : 3.0) * theta;
  const auto nph   = box.dE(ell) / (xmean * wph);
  auto rng0        = runko::Rng { runko::rng_key(par.seed, istep_, 0, stream::photons_count), 0 };
  const auto n     = static_cast<std::size_t>(std::floor(nph + rng0.uniform()));
  if(n == 0) { return; }

  Ctx c { t, par.nbin };
  const auto stage = this->stager(n);
  tally_[0] = Tally {};

  const auto db   = ptr(tally_);
  const auto key  = runko::rng_key(par.seed, istep_, 0, stream::photons);
  const auto mins = t.mins, maxs = t.maxs;
  const auto lay = c.v[0].lay;
  const auto w   = tyvi::mdgrid_work {};
  const auto i0  = thrust::counting_iterator<I>(0);

  // one lane per photon: uniform position, energy from the spectrum, isotropic
  thrust::for_each(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    [=](const I i) {
      auto rng           = runko::Rng { key, i };
      const float pos[3] = {
        mins[0] + (maxs[0] - mins[0]) * static_cast<float>(rng.uniform()),
        mins[1] + (maxs[1] - mins[1]) * static_cast<float>(rng.uniform()),
        mins[2] + (maxs[2] - mins[2]) * static_cast<float>(rng.uniform())
      };

      // Wien: a Gamma(3, Theta) deviate; Planck: first the term j of the Bose series
      // (weight 1/j^3, total zeta(3)), then Gamma(3, Theta)/j
      auto j = 1.0;
      if(planck) {
        const auto z1 = 1.2020569 * rng.uniform();
        auto cum      = 1.0;
        while(z1 > cum and j < 50.0) {
          j += 1.0;
          cum += 1.0 / (j * j * j);
        }
      }
      const auto x =
        -theta * sstd::log(rng.uniform() * rng.uniform() * rng.uniform()) / j;

      // drawn in this order: direction, then sub-range
      const auto k   = x * iso_dir(rng);
      const auto bit = static_cast<int>(rng.uniform() < 0.5);
      stage(k, pos, wph, 1, lay.cell_of(pos[0], pos[1], pos[2]), Species::ph, bit);
      sstd::atomic_add(&db->E, wph * x);
    });
  w.wait();

  {
    const Tally h = tally_[0];
    budget.inj_E += h.E;
  }

  this->place(c);
}

void
  Sampler::heat_diffusive(TileView& t, const double ell)
{
  this->lay_if_needed(t);

  // Ito kick of the isotropic Fokker-Planck operator (1/p^2) d/dp [p^2 D dn/dp], taken
  // on the momentum vector so that it also scatters the direction. Expected gain per
  // lepton: D_0 p^2 (4 gamma^2 + 1)/gamma^3; D_0 is set every step so that the summed
  // gain is the energy the compactness injects.
  Ctx c { t, par.nbin };
  const auto w = tyvi::mdgrid_work {};

  // tot = sum of w p^2 (4 gamma^2 + 1)/gamma^3 over the live leptons
  auto tot = 0.0;
  for(auto s = 1uz; s <= 2uz; ++s) {
    const auto v  = c.v[s];
    const auto n  = t.con[s]->size();
    const auto i0 = thrust::counting_iterator<I>(0);
    tot += thrust::transform_reduce(
      w.on_this(),
      i0,
      i0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return 0.0; }
        const auto p   = v.u(i);
        const auto p2  = toolbox::dot(p, p);
        const auto gam = sstd::sqrt(1.0 + p2);
        return v.w(i) * p2 * (4.0 * gam * gam + 1.0) / (gam * gam * gam);
      },
      0.0,
      thrust::plus<double> {});
  }

  if(tot <= 0.0) { return; }
  const auto D0 = box.dE(ell) / tot;
  tally_[0] = Tally {};
  const auto db = ptr(tally_);

  // the kicks: one lane per lepton slot of either species
  for(auto s = 1uz; s <= 2uz; ++s) {
    const auto v   = c.v[s];
    const auto n   = t.con[s]->size();
    const auto key = runko::rng_key(par.seed, istep_, s, stream::heating);
    const auto i0  = thrust::counting_iterator<I>(0);

    thrust::for_each(
      w.on_this(),
      i0,
      i0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return; }

        auto rng      = runko::Rng { key, i };
        const auto p  = v.u(i);
        const auto p2 = toolbox::dot(p, p);
        const auto sd = sstd::sqrt(2.0 * D0 * p2);  // width of the Ito kick

        const auto p1 = (1.0 + 2.0 * D0) * p +
                        Vec3(sd * rng.normal(), sd * rng.normal(), sd * rng.normal());

        sstd::atomic_add(
          &db->E,
          v.w(i) * (sstd::sqrt(1.0 + toolbox::dot(p1, p1)) - sstd::sqrt(1.0 + p2)));
        v.set_u(i, p1);
        v.touch(i);
      });
  }
  w.wait();

  const Tally h = tally_[0];
  budget.heat_E += h.E;
}

void
  Sampler::escape_photons(TileView& t, const double tau_ext, const double dt_ext)
{
  this->lay_if_needed(t);
  Ctx c { t, par.nbin };

  // pair depth and step in slab crossings: the slab's, given by the driver, or this
  // tile's own
  const auto tau_T  = tau_ext >= 0.0
                        ? tau_ext
                        : t.layout[1]->totals(*t.con[1]).W + t.layout[2]->totals(*t.con[2]).W;
  const auto dt_esc = dt_ext > 0.0 ? dt_ext : box.dt;
  if(esc_hist_.size() != static_cast<std::size_t>(n_esc)) {
    esc_hist_.assign(n_esc, 0.0);
    esc_mu_hist_.assign(static_cast<std::size_t>(n_esc) * n_mu, 0.0);
  }

  tally_[0] = Tally {};
  const auto db   = ptr(tally_);
  const auto hist = thrust::raw_pointer_cast(esc_hist_.data());
  const auto hmu  = thrust::raw_pointer_cast(esc_mu_hist_.data());

  const auto v   = c.v[0];
  const auto n   = t.con[0]->size();
  const auto key = runko::rng_key(par.seed, istep_, 0, stream::escape);
  const auto llo = std::log(esc_lo), dl = (std::log(esc_hi) - std::log(esc_lo)) / n_esc;
  const auto w  = tyvi::mdgrid_work {};
  const auto i0 = thrust::counting_iterator<I>(0);

  // one lane per photon slot: escape with the slab probability, book x ell_x per log
  // bin and per (log bin, cosine to z)
  thrust::for_each(
    w.on_this(),
    i0,
    i0 + static_cast<std::ptrdiff_t>(n),
    [=](const I i) {
      if(not v.alive(i)) { return; }
      auto rng     = runko::Rng { key, i };
      const auto x = v.e(i);
      if(not(rng.uniform() < Box::p_esc(x, tau_T, dt_esc))) { return; }

      const auto wg = v.w(i);
      v.kill(i);
      sstd::atomic_add(&db->E, x * wg);
      sstd::atomic_add(&db->W, wg);

      const auto b = sstd::clamp(
        static_cast<int>(sstd::floor((sstd::log(x) - llo) / dl)),
        0,
        n_esc - 1);
      sstd::atomic_add(&hist[b], x * wg);

      const auto mu = static_cast<double>(v.vel[i][2]) / x;
      const auto m =
        sstd::clamp(static_cast<int>(sstd::floor((mu + 1.0) * 0.5 * n_mu)), 0, n_mu - 1);
      sstd::atomic_add(&hmu[b * n_mu + m], x * wg);
    });
  w.wait();

  const Tally h = tally_[0];
  budget.esc_E += h.E;
  budget.esc_W += h.W;
}

// zeros before the first escape
std::vector<double>
  Sampler::escaped_mu_spectrum() const
{
  if(esc_mu_hist_.empty()) {
    return std::vector<double>(static_cast<std::size_t>(n_esc) * n_mu, 0.0);
  }
  thrust::host_vector<double> h = esc_mu_hist_;
  return std::vector<double>(h.begin(), h.end());
}

// zeros before the first escape
std::vector<double>
  Sampler::escaped_spectrum() const
{
  if(esc_hist_.empty()) { return std::vector<double>(n_esc, 0.0); }
  thrust::host_vector<double> h = esc_hist_;
  return std::vector<double>(h.begin(), h.end());
}

// start a new accounting window: budgets and the escaping spectrum
void
  Sampler::reset_budget()
{
  budget = Budget {};

  const auto w = tyvi::mdgrid_work {};
  thrust::fill(w.on_this(), esc_hist_.begin(), esc_hist_.end(), 0.0);
  thrust::fill(w.on_this(), esc_mu_hist_.begin(), esc_mu_hist_.end(), 0.0);
  w.wait();
}

double
  Sampler::tau_ssa(const TileView& t, const double x) const
{
  // With weight in units of depth, a photon's absorptions per crossing are sum w 2
  // sigma-hat F over the leptons: the box depth in this channel, the analogue of tau_T.
  // Its unit crossing marks x_SSA.
  if(box.b_perp() <= 0.0) { return 0.0; }

  // each lepton absorbs at the tile field across its pitch angle
  const auto fv = t.field;
  const auto w  = tyvi::mdgrid_work {};
  auto out      = 0.0;
  const auto k  = Vec3(x, 0.0, 0.0);  // the rate depends on x, gamma and the field only

  for(auto s = 1uz; s <= 2uz; ++s) {
    auto& con     = *t.con[s];
    const auto v  = make_view(con, box.w_lep);
    const auto n  = con.size();
    const auto i0 = thrust::counting_iterator<I>(0);
    out += thrust::transform_reduce(
      w.on_this(),
      i0,
      i0 + static_cast<std::ptrdiff_t>(n),
      [=](const I i) {
        if(not v.alive(i)) { return 0.0; }

        const auto u = v.u(i);
        const auto b = fv.b_ssa(u / sstd::max(toolbox::norm(u), 1e-300));
        return v.w(i) * 2.0 * ssa::rate(k, u, b, 0.0);
      },
      0.0,
      thrust::plus<double> {});
  }

  return out;
}

}  // namespace qed
