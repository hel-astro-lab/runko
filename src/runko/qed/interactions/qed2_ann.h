// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Pair annihilation, e- + e+ -> photon + photon: the Dirac cross section through its
// energy bracket, the rate (sigma_ann times the Moller relative velocity, halved), its
// bound over a pair of energy bins, and the product sampler. Units: m_e c^2, m_e c.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::ann {

inline constexpr double SHMAX = 0.208;  // sup of rate, 0.20798 at gamma_cm = 1.1354

/// Returns the energy bracket of the pair invariant q_e = p+ . p-, with
/// sigma_ann F/2 = (3/32) bracket(q_e)/(g- g+); increasing in q_e.
[[nodiscard]]
constexpr double
  bracket(double qe)
{
  qe            = sstd::max(qe, 1.0 + 1e-14);  // q_e = 1: co-moving pair
  const auto g2 = (qe + 1.0) / 2.0;            // gamma_cm^2
  const auto b  = sstd::sqrt(1.0 - 1.0 / g2);

  // 2 artanh(b)/b, with its slow-pair series 2 + 2b^2/3
  const auto Lb = b > 1e-3 ? 2.0 * sstd::atanh(sstd::min(b, 1.0 - 1e-16)) / b
                           : 2.0 + 2.0 * b * b / 3.0;
  return (2.0 + 2.0 / g2 - 1.0 / (g2 * g2)) * Lb - 2.0 - 2.0 / g2;
}

/// Returns sigma_ann F_Moller/2 = (3/32)/(g+ g-) bracket(q_e) for e- um, e+ up.
[[nodiscard]]
constexpr double
  rate(const Vec3& um, const Vec3& up)
{
  // the Moller factor cancels the 1/(beta_cm gamma_cm^2) of the Dirac cross section
  const auto gm = sstd::sqrt(1.0 + toolbox::dot(um, um));
  const auto gp = sstd::sqrt(1.0 + toolbox::dot(up, up));
  const auto qe = gm * gp - toolbox::dot(um, up);
  return (3.0 / 32.0) / (gm * gp) * bracket(qe);
}

/// Returns a bound on rate over e- in [g1lo, g1hi] and e+ in [g2lo, g2hi].
[[nodiscard]]
constexpr double
  bound(const double g1lo, const double g1hi, const double g2lo, const double g2hi)
{
  // largest q_e (head-on at the top corner) over the smallest g- g+
  const auto q  = g1hi * g2hi + sstd::sqrt(sstd::max(g1hi * g1hi - 1.0, 0.0)) *
                                  sstd::sqrt(sstd::max(g2hi * g2hi - 1.0, 0.0));
  const auto sh = (3.0 / 32.0) * bracket(q) / (g1lo * g2lo);
  return sh == sh and sh < 1e300 ? sstd::min(sh, SHMAX) : 0.0;
}

struct [[nodiscard]] Products {
  Vec3 k0, k1;
};

/// Returns the photons of e- um, e+ up; ntry receives the angular-draw proposals.
[[nodiscard]]
constexpr Products
  sample(const Vec3& um, const Vec3& up, runko::Rng& rng, int& ntry)
{
  // the photons leave the CM back to back with energy gamma_cm each
  const auto gm   = sstd::sqrt(1.0 + toolbox::dot(um, um));
  const auto gp   = sstd::sqrt(1.0 + toolbox::dot(up, up));
  const auto s0   = gm + gp;
  const auto svec = um + up;
  const auto q    = gp * gm - toolbox::dot(up, um) + 1.0;  // q_e + 1 = 2 gamma_cm^2
  const auto bc   = svec / (-s0);                          // minus the CM velocity
  const auto vc   = toolbox::norm(bc);
  const auto gcm  = sstd::sqrt(q / 2.0);
  const auto vcm  = sstd::sqrt(sstd::max(1.0 - 1.0 / (gcm * gcm), 0.0));
  const auto nm   = toolbox::norm(um);
  const auto omm  = nm > 1e-300 ? um / nm : Vec3(1.0, 0.0, 0.0);

  // polar axis: the boost, or the electron direction (y = 1) for a pair at rest;
  // y = cosine of the electron's CM direction to the axis
  const auto deg  = vc < 1e-10;
  const auto kvec = deg ? omm : bc / sstd::max(vc, 1e-300);
  auto y = deg
             ? 1.0
             : (gp - gm) / ((gp + gm) * sstd::max(vc, 1e-300) * sstd::max(vcm, 1e-300));
  y      = sstd::clamp(y, -1.0, 1.0);

  const auto t   = perp_basis(kvec, omm);
  const auto omr = cm_direction(y, gcm, vcm, t, rng, ntry);

  const auto A = boost(bc, gcm, gcm * omr);
  const auto B = boost(bc, gcm, (-gcm) * omr);
  return Products { A.p, B.p };
}

}  // namespace qed::ann
