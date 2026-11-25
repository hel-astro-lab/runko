// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/qed/interactions/qed1_syn.h"
#include "runko/tools/math.h"
#include "runko/tools/random.h"

namespace qed {

/// Exchanges a light partner of energy gam may take as one.
///
/// It is kappa (gam - 1) / dgmax, dgmax being the most energy a single exchange gives
/// it; at least 1, and 1 for kappa <= 0.
[[nodiscard]]
constexpr double
  fkin(const double kappa, const double gam, const double dgmax)
{
  if(kappa <= 0.0) { return 1.0; }

  return sstd::max(1.0, kappa * (gam - 1.0) / sstd::max(dgmax, 1e-300));
}

/// Weight S an event of the pair moves.
///
/// Two LPs of weights w_i and w_j stand for w_i w_j 2 r / N_2 reactions per step, r
/// being the rate function of the channel. An event moves the weight S and fires with
/// the probability (w_i w_j / S) 2 r / N_2, so that probability times size is the
/// reaction count whatever S is chosen. The incident is updated with the probability
/// S / w_i; the partner takes S / w_j exchanges, a probability while below one and a
/// multiplicity above. S is min(w_i, w_j) between equals (fk = 1), and up to the whole
/// of a heavy incident when the light partner takes fk exchanges as one.
[[nodiscard]]
constexpr double
  event_size(const double wi, const double wj, const double fk)
{ return sstd::min(wi, fk * wj); }

/// Weight of the pair in its event probability, w_i w_j / event_size.
///
/// Rising in both weights and falling in fk, so the largest weights of two groups and
/// their smallest fk bound it.
[[nodiscard]]
constexpr double
  pair_weight(const double wi, const double wj, const double fk)
{ return sstd::max(wj, wi / fk); }


/// Photon resolution function: LPs per unit weight go as (x/x_ref)^a.
[[nodiscard]]
constexpr double
  f_res(const double e, const double a, const double x_ref)
{ return sstd::pow(e / x_ref, a); }

/// A product realized as n_r copies of weight w_each.
struct Copies {
  int n_r;
  double w_each;
};

/// Realize the copy number n of a product of weight s as n_r copies of w_each.
///
/// n >= 1 splits into floor(n) or floor(n) + 1 copies of s/n_r; n < 1 keeps one copy
/// of s/n with probability n (Russian roulette). Both keep <n_r> = min(n, n_cap).
[[nodiscard]]
constexpr Copies
  copies(double n, const double s, runko::Rng& rng, const int n_cap)
{
  n = sstd::min(n, static_cast<double>(n_cap));

  const auto n_r = static_cast<int>(sstd::floor(n + rng.uniform()));
  const auto den =
    n >= 1.0 ? static_cast<double>(sstd::max(n_r, 1)) : sstd::max(n, 1e-300);
  return Copies { n_r, s / den };
}


/// Accumulation factor of a resolved emission of mean photon energy xm.
///
/// It is acc_syn max(1, xm^-2), capped at kappa (gamma - 1)/xm, at least 1.
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

/// Parent energy after one step of synchrotron cooling above the floor.
///
/// The cooling dgamma/dt = -(2/3)(alpha c/lamC) chi^2 g Phi is integrated in closed
/// form for chi ~ gamma. Phi is the power fraction above the floor, from its rate lam
/// and mean energy xm.
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

/// What a lepton of energy gam does in a step in which it emits lam photons of mean
/// energy xm, before any photon is drawn.
struct EmissionPlan {
  /// Below one emission per step.
  bool resolved;
  /// An LP is emitted this step.
  bool fire;
  /// Emissions merged into the LP; 1 if unresolved.
  double f;
  /// Parent energy after the continuous cooling; gam if resolved.
  double g_cooled;
};

/// The emission plan of a lepton for one step; z is the uniform that decides.
///
/// Resolved (lam <= 1): one LP stands for f emissions and fires with the probability
/// lam / f <= 1, so that the emissions of a step average lam. Unresolved: the parent
/// cools along gamma_cooled whatever is emitted, and one LP stands for all the
/// emissions of the step.
[[nodiscard]]
constexpr EmissionPlan
  emission_plan(
    const double gam,
    const double lam,
    const double xm,
    const double b,
    const double N1,
    const double acc_syn,
    const double kappa,
    const double z)
{
  const auto resolved = lam <= 1.0;
  const auto f        = resolved ? facc_syn(gam, xm, acc_syn, kappa) : 1.0;
  return EmissionPlan { .resolved = resolved,
                        .fire     = resolved ? z < lam / f : true,
                        .f        = f,
                        .g_cooled =
                          resolved ? gam : sstd::max(gamma_cooled(gam, lam, xm, b, N1), 1.0) };
}

/// The emitted LP and its parent once the photon energy x is drawn.
struct Emitted {
  /// Parent energy after the step.
  double g1;
  /// Weight of the photon LP.
  double wem;
};

/// The emitted LP of plan p and its parent of weight w once x is drawn.
///
/// A resolved parent recoils by what its LP carries, f x; an unresolved one has
/// cooled, and its LP stands for loss / xm emissions, so that the emitted energy is
/// the loss on average. A recoil is clipped at the mass shell and the emitted weight
/// shrinks with it.
[[nodiscard]]
constexpr Emitted
  emitted(
    const EmissionPlan& p,
    const double w,
    const double gam,
    const double xm,
    const double x)
{
  const auto dgu  = gam - p.g_cooled;  // energy the continuous cooling took
  const auto mult = p.resolved ? 1.0 : dgu / sstd::max(xm, 1e-300);
  const auto dg   = p.resolved ? p.f * x : dgu;
  const auto dgc  = sstd::min(dg, gam - 1.0);
  const auto wem  = w * p.f * mult * (dg > 0.0 ? dgc / sstd::max(dg, 1e-300) : 1.0);
  return Emitted { .g1 = p.resolved ? sstd::max(gam - dgc, 1.0) : p.g_cooled, .wem = wem };
}

}  // namespace qed
