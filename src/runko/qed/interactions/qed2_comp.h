// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/qed/interactions/qed_kinematics.h"

/// Compton scattering, photon + lepton -> photon + lepton.
///
/// The rate is the Klein-Nishina cross section times the Moller relative velocity,
/// halved. Units: m_e c^2, m_e c, sigma_T.
namespace qed::comp {

/// Supremum of rate: sigma_KN -> 1 as x_R -> 0.
inline constexpr double rate_max = 1.0;

/// Returns the Klein-Nishina cross section from the rest-frame photon energy x_R.
[[nodiscard]]
constexpr double
  sigma_kn(double xR)
{
  // The closed form cancels to 0/0 as x_R -> 0: series below 1e-2
  xR = sstd::max(xR, 1e-12);
  if(xR < 1e-2) { return 1.0 - 2.0 * xR + 5.2 * xR * xR - 13.3 * xR * xR * xR; }

  const auto l = sstd::log1p(2.0 * xR);
  return 0.375 / (xR * xR) *
         (4.0 + (xR - 2.0 - 2.0 / xR) * l +
          2.0 * xR * xR * (1.0 + xR) / ((1.0 + 2.0 * xR) * (1.0 + 2.0 * xR)));
}

/// Returns sigma_KN(x_R) (1 - beta mu)/2 for photon k, lepton u; x_R = gamma x - u.k.
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const Vec3& u)
{
  const auto x  = toolbox::norm(k);
  const auto g  = sstd::sqrt(1.0 + toolbox::dot(u, u));
  const auto xR = sstd::max(g * x - toolbox::dot(u, k), 1e-12);
  const auto F  = xR / (g * x);  // 1 - beta mu
  return 0.5 * sigma_kn(xR) * F;
}

/// Returns a bound on rate over photons in [xlo, xhi] and leptons in [glo, ghi], 0
/// for an empty bin.
[[nodiscard]]
constexpr double
  bound(const double xlo, const double xhi, const double glo, const double ghi)
{
  // The rate sigma_KN x_R/2 grows with x_R: take the largest x_R (head-on at the top
  // corner) over the smallest x gamma.
  const auto y  = xhi * (ghi + sstd::sqrt(sstd::max(ghi * ghi - 1.0, 0.0)));
  const auto sh = sigma_kn(y) * y / 2.0 / (xlo * glo);
  return sh == sh and sh < 1e300 ? sstd::min(sh, rate_max) : 0.0;  // NaN/inf: empty bin
}

/// Returns the largest energy change of the lepton per event, 4 x gamma^2.
[[nodiscard]]
constexpr double
  dgmax(const double x, const double g)
{ return 4.0 * x * g * g; }

struct Products {
  /// Scattered photon.
  Vec3 k1;
  /// Recoiled lepton.
  Vec3 u1;
};

/// Klein-Nishina scattering of photon k on lepton u; ntry receives the proposals used.
[[nodiscard]]
constexpr Products
  sample(const Vec3& k, const Vec3& u, runko::Rng& rng, int& ntry)
{
  // Scatter in the lepton rest frame, boost back, recoil from momentum conservation.
  // Draw t = x_R/x'_R on [1, 1 + 2 x_R] from f(t) = 1/t^3 + 1/t + (mu^2 - 1)/t^2,
  // mu = 1 - (t-1)/x_R, bounded by 2/t: propose from 1/t, accept with t f(t)/2;
  // acceptance 0.5 to 0.67 at every incident energy.
  const auto x0  = toolbox::norm(k);
  const auto g0  = sstd::sqrt(1.0 + toolbox::dot(u, u));
  const auto b0  = u / g0;
  const auto R   = boost(b0, x0, k);
  const auto xR  = R.e;
  const auto omR = R.p / xR;

  // phi is uniform: any azimuth origin not parallel to omR serves
  const auto nb  = toolbox::norm(b0);
  const auto ref = nb > 1e-12 ? b0 / nb : Vec3(1.0, 0.0, 0.0);
  const auto t   = perp_basis(omR, ref);

  const auto a   = sstd::max(xR, 1e-300);
  const auto lnT = sstd::log1p(2.0 * a);
  auto mu        = 1.0;
  ntry           = 0;

  while(ntry < trip_cap) {
    ++ntry;
    const auto tm1 = sstd::expm1(rng.uniform() * lnT);  // t - 1, exact for soft photons
    const auto tt  = 1.0 + tm1;
    mu             = 1.0 - tm1 / a;
    if((1.0 / (tt * tt) + 1.0 + (mu * mu - 1.0) / tt) / 2.0 > rng.uniform()) { break; }
  }

  const auto x1R = xR / (1.0 + xR * (1.0 - mu));
  const auto phi = sstd::TWOPI * rng.uniform();
  const auto sz  = sstd::sqrt(sstd::max(1.0 - mu * mu, 0.0));
  const auto om1R =
    (sz * sstd::sin(phi)) * t.i + (sz * sstd::cos(phi)) * t.j + mu * t.k;

  const auto L = boost((-1.0) * b0, x1R, x1R * om1R);
  return Products { L.p, u + k - L.p };
}

}  // namespace qed::comp
