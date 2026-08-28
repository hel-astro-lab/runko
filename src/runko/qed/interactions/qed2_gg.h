// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Photon-photon pair production, photon + photon -> e- + e+: the Breit-Wheeler cross
// section, the reaction rate built on it, its bound over a pair of energy bins, and the
// product sampler. The threshold is imposed by the cross section itself, which is exactly
// zero below it. Units: m_e c^2, m_e c, sigma_T.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::gg {

inline constexpr double SHMAX = 0.256;  // bounds rate, whose supremum is 0.25564 at q = 3.94

/// Breit-Wheeler total cross section in Thomson units, from q = x1 x2 (1 - mu) = 2 gamma_cm^2.
[[nodiscard]]
constexpr double
  sigma_gg(const double q)
{
  if(not(q > 2.0 * (1.0 + 1e-12))) { return 0.0; }         // below the pair threshold
  const auto b  = sstd::sqrt(1.0 - 2.0 / q);                // beta_cm of each lepton
  const auto Lb = 2.0 * sstd::atanh(sstd::min(b, 1.0 - 1e-16));  // ln[(1+b)/(1-b)]
  return (3.0 / 16.0) * (1.0 - b * b) * ((3.0 - b * b * b * b) * Lb - 2.0 * b * (2.0 - b * b));
}

/// photons k1, k2 -> sigma_gg (1 - mu)/2; zero below threshold x1 x2 (1 - mu) < 2
[[nodiscard]]
constexpr double
  rate(const Vec3& k1, const Vec3& k2)
{
  const auto x1 = toolbox::norm(k1);
  const auto x2 = toolbox::norm(k2);
  const auto q  = x1 * x2 - toolbox::dot(k1, k2);           // x1 x2 (1 - mu)
  return 0.5 * sigma_gg(q) * q / sstd::max(x1 * x2, 1e-300);  // F = 1 - mu = q/(x1 x2)
}

/// Bound on rate over photons in [x1lo, x1hi] and [x2lo, x2hi]; exactly 0 below threshold.
[[nodiscard]]
constexpr double
  bound(const double x1lo, const double x1hi, const double x2lo, const double x2hi)
{
  const auto q  = 2.0 * x1hi * x2hi;                        // largest invariant, head-on
  const auto sh = sigma_gg(q) * q / 2.0 / (x1lo * x2lo);
  return sh == sh and sh < 1e300 ? sstd::min(sh, SHMAX) : 0.0;
}

struct [[nodiscard]] Products {
  Vec3 p1, p2;  // the two leptons; the caller assigns the charges
};

/// photons k1, k2 -> lepton momenta; ntry receives the proposals of the angular draw.
[[nodiscard]]
constexpr Products
  sample(const Vec3& k1, const Vec3& k2, runko::Rng& rng, int& ntry)
{
  // Crossing-symmetric counterpart of the annihilation sampler: in the CM the two leptons
  // leave back to back with energy gamma_cm and momentum gamma_cm beta_cm.
  const auto x1   = toolbox::norm(k1);
  const auto x2   = toolbox::norm(k2);
  const auto s0   = x1 + x2;                                // total energy
  const auto svec = k1 + k2;                                // total momentum
  const auto s    = toolbox::norm(svec);
  const auto q    = x1 * x2 - toolbox::dot(k1, k2);         // 2 gamma_cm^2
  const auto bc   = svec / (-s0);                           // minus the CM velocity
  const auto gcm  = sstd::sqrt(q / 2.0);
  const auto vcm  = sstd::sqrt(sstd::max(1.0 - 1.0 / (gcm * gcm), 0.0));
  const auto om1  = k1 / x1;                                // direction of photon 1
  const auto deg  = s < 1e-10;                              // head-on and equal: no boost axis
  const auto kvec = deg ? om1 : bc / sstd::max(toolbox::norm(bc), 1e-300);
  const auto cz1  = deg ? 1.0 : sstd::clamp((x2 - x1) / sstd::max(s, 1e-300), -1.0, 1.0);
  const auto t    = perp_basis(kvec, om1);
  const auto omr  = cm_direction(cz1, gcm, vcm, t, rng, ntry);  // lepton direction in the CM
  const auto A    = boost(bc, gcm, (gcm * vcm) * omr);
  const auto B    = boost(bc, gcm, (-gcm * vcm) * omr);
  return Products { A.p, B.p };
}

}  // namespace qed::gg
