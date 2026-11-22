// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Compton scattering, photon + lepton -> photon + lepton: the Klein-Nishina cross
// section, the reaction rate built on it, its bound over a pair of energy bins, and the
// product sampler. rate() is sigma_KN times the Moller relative velocity, halved, so
// that it is directly comparable with the constant bound SHMAX. Units: m_e c^2, m_e c,
// sigma_T.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::comp {

inline constexpr double SHMAX = 1.0;  // sup of rate: sigma_KN -> 1 as x_R -> 0

/// Klein-Nishina total cross section in Thomson units, from the rest-frame photon
/// energy.
[[nodiscard]]
constexpr double
  sigma_kn(double xR)
{
  // Tends to 1 in the Thomson limit and falls as 1/x_R once the recoil matters. The
  // closed form cancels to 0/0 as x_R -> 0, so small arguments take its series instead.
  xR = sstd::max(xR, 1e-12);  // guard the 1/x_R^2 of the closed form
  if(xR < 1e-2) {
    return 1.0 - 2.0 * xR + 5.2 * xR * xR - 13.3 * xR * xR * xR;
  }                                      // series
  const auto l = sstd::log1p(2.0 * xR);  // ln(1 + 2 x_R) of the closed form
  return 0.375 / (xR * xR) *
         (4.0 + (xR - 2.0 - 2.0 / xR) * l +
          2.0 * xR * xR * (1.0 + xR) / ((1.0 + 2.0 * xR) * (1.0 + 2.0 * xR)));
}

/// photon k, lepton u -> sigma_KN(x_R) (1 - beta mu)/2, with x_R = gamma x - u.k
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const Vec3& u)
{
  const auto x = toolbox::norm(k);                      // photon energy
  const auto g = sstd::sqrt(1.0 + toolbox::dot(u, u));  // lepton Lorentz factor
  const auto xR =
    sstd::max(g * x - toolbox::dot(u, k), 1e-12);  // = gamma x (1 - beta mu)
  const auto F = xR / (g * x);                     // Moller factor 1 - beta mu
  return 0.5 * sigma_kn(xR) * F;                   // bounded by SHMAX
}

/// Bound on rate over photons in [xlo, xhi] and leptons in [glo, ghi]; 0 for an empty
/// bin.
[[nodiscard]]
constexpr double
  bound(const double xlo, const double xhi, const double glo, const double ghi)
{
  // The rate is monotone in the rest-frame energy over a product of the two energies,
  // so its bound over a bin pair is the largest rest-frame energy the bins can reach
  // over the smallest energy product: no search, no table. largest x_R: head-on at the
  // top corner
  const auto y = xhi * (ghi + sstd::sqrt(sstd::max(ghi * ghi - 1.0, 0.0)));
  const auto sh =
    sigma_kn(y) * y / 2.0 / (xlo * glo);  // sigma_KN x_R/2 over the smallest x gamma
  return sh == sh and sh < 1e300 ? sstd::min(sh, SHMAX) : 0.0;  // NaN/inf: empty bin
}

/// Largest energy change of the lepton per event, for the accumulation factor.
[[nodiscard]]
constexpr double
  dgmax(const double x, const double g)
{ return 4.0 * x * g * g; }

struct [[nodiscard]] Products {
  Vec3 k1;  // scattered photon
  Vec3 u1;  // recoiled lepton
};

/// Klein-Nishina scattering of photon k on lepton u; ntry receives the proposals used.
[[nodiscard]]
constexpr Products
  sample(const Vec3& k, const Vec3& u, runko::Rng& rng, int& ntry)
{
  // The angular distribution is simple only in the rest frame of the lepton, so the
  // photon is boosted there, scattered, and boosted back; the lepton follows from
  // four-momentum conservation. What is drawn is the energy ratio t = x_R/x'_R on [1, 1
  // + 2 x_R], density
  //     f(t) = 1/t^3 + 1/t + (mu^2 - 1)/t^2,   mu = 1 - (t-1)/x_R,
  // bounded by 2/t. Proposing t from 1/t and thinning with t f(t)/2 accepts 0.52 to
  // 0.67 of the proposals at every incident energy.
  const auto x0 = toolbox::norm(k);  // incoming photon energy
  const auto g0 =
    sstd::sqrt(1.0 + toolbox::dot(u, u));  // incoming lepton Lorentz factor
  const auto b0  = u / g0;                 // its velocity: boost to rest frame
  const auto R   = boost(b0, x0, k);       // photon in the lepton rest frame
  const auto xR  = R.e;
  const auto omR = R.p / xR;  // its direction there
  // Azimuth origin: the plane of the photon direction and the lepton velocity; any
  // plane does, phi is uniform, but ref must not be parallel to omR.
  const auto nb  = toolbox::norm(b0);
  const auto ref = nb > 1e-12 ? b0 / nb : Vec3(1.0, 0.0, 0.0);
  const auto t   = perp_basis(omR, ref);  // incoming photon as polar axis
  const auto a   = sstd::max(xR, 1e-300);
  const auto lnT = sstd::log1p(2.0 * a);  // log of the largest ratio 1 + 2 x_R
  auto mu        = 1.0;                   // cos theta in the rest frame
  ntry           = 0;
  // rejection loop: t from 1/t on [1, 1 + 2 x_R], accept with t f(t)/2
  while(ntry < TRIP_CAP) {
    ++ntry;
    const auto tm1 =
      sstd::expm1(rng.uniform() * lnT);  // t - 1, exact for a soft photon
    const auto tt = 1.0 + tm1;
    mu            = 1.0 - tm1 / a;  // cos theta in the rest frame
    if((1.0 / (tt * tt) + 1.0 + (mu * mu - 1.0) / tt) / 2.0 > rng.uniform()) { break; }
  }
  const auto x1R = xR / (1.0 + xR * (1.0 - mu));  // scattered energy in the rest frame
  const auto phi = sstd::TWOPI * rng.uniform();
  const auto sz  = sstd::sqrt(sstd::max(1.0 - mu * mu, 0.0));  // sin(theta_R)
  // scattered direction, rest frame
  const auto om1R =
    (sz * sstd::sin(phi)) * t.i + (sz * sstd::cos(phi)) * t.j + mu * t.k;
  const auto L = boost((-1.0) * b0, x1R, x1R * om1R);  // back in the laboratory
  return Products { L.p, u + k - L.p };  // lepton recoil from conservation
}

}  // namespace qed::comp
