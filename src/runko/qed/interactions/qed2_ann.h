// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Pair annihilation, e- + e+ -> photon + photon: the Dirac cross section in the form
// its energy bracket takes, the reaction rate built on it, its bound over a pair of
// energy bins, and the product sampler. rate() is sigma_ann times the Moller relative
// velocity, halved, comparable with the constant bound SHMAX. Units: m_e c^2, m_e c,
// sigma_T.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::ann {

inline constexpr double SHMAX = 0.208;  // sup of rate, 0.20798 at gamma_cm = 1.1354

/// Energy bracket of the annihilation rate: sigma_ann F/2 = (3/32) bracket(q_e)/(g-
/// g+).
[[nodiscard]]
constexpr double
  bracket(double qe)
{
  // A function of the pair invariant q_e = p_+ . p_- alone, increasing in it like ln
  // q_e, which is what lets a bound on q_e alone bound the rate of a whole bin pair.
  qe = sstd::max(qe, 1.0 + 1e-14);   // q_e = 1 is the degenerate co-moving pair
  const auto g2 = (qe + 1.0) / 2.0;  // gamma_cm^2
  const auto b  = sstd::sqrt(1.0 - 1.0 / g2);  // beta_cm
  // L(b)/b = 2 artanh(b)/b, with its slow-pair series 2 + 2b^2/3
  const auto Lb = b > 1e-3 ? 2.0 * sstd::atanh(sstd::min(b, 1.0 - 1e-16)) / b
                           : 2.0 + 2.0 * b * b / 3.0;
  return (2.0 + 2.0 / g2 - 1.0 / (g2 * g2)) * Lb - 2.0 - 2.0 / g2;
}

/// e- um, e+ up -> sigma_ann F_Moller/2 = (3/32)/(g+ g-) bracket(q_e)
[[nodiscard]]
constexpr double
  rate(const Vec3& um, const Vec3& up)
{
  // The Moller factor of a lepton pair cancels the 1/(beta_cm gamma_cm^2) prefactor of
  // the Dirac cross section exactly, so only 3/32 over g+ g- survives in front of the
  // bracket.
  const auto gm = sstd::sqrt(1.0 + toolbox::dot(um, um));  // electron Lorentz factor
  const auto gp = sstd::sqrt(1.0 + toolbox::dot(up, up));  // positron Lorentz factor
  const auto qe = gm * gp - toolbox::dot(um, up);          // pair invariant q_e >= 1
  return (3.0 / 32.0) / (gm * gp) * bracket(qe);
}

/// Bound on rate over electrons in [g1lo, g1hi] and positrons in [g2lo, g2hi].
[[nodiscard]]
constexpr double
  bound(const double g1lo, const double g1hi, const double g2lo, const double g2hi)
{
  const auto q  = g1hi * g2hi + sstd::sqrt(sstd::max(g1hi * g1hi - 1.0, 0.0)) *
                                  // largest q_e: head-on at the top corner
                                  sstd::sqrt(sstd::max(g2hi * g2hi - 1.0, 0.0));
  const auto sh = (3.0 / 32.0) * bracket(q) / (g1lo * g2lo);  // over the smallest g- g+
  return sh == sh and sh < 1e300 ? sstd::min(sh, SHMAX) : 0.0;
}

struct [[nodiscard]] Products {
  Vec3 k0, k1;  // the two photons
};

/// e- um, e+ up -> two photon momenta; ntry receives the proposals of the angular draw.
[[nodiscard]]
constexpr Products
  sample(const Vec3& um, const Vec3& up, runko::Rng& rng, int& ntry)
{
  // In the CM frame the two photons leave back to back with energy gamma_cm each, so
  // only one direction is drawn; the pair kinematics fix everything else.
  // ---- the pair in the laboratory and its CM frame
  const auto gm   = sstd::sqrt(1.0 + toolbox::dot(um, um));  // electron Lorentz factor
  const auto gp   = sstd::sqrt(1.0 + toolbox::dot(up, up));  // positron Lorentz factor
  const auto s0   = gm + gp;                                 // total energy
  const auto svec = um + up;                                 // total momentum
  const auto q    = gp * gm - toolbox::dot(up, um) + 1.0;    // q_e + 1 = 2 gamma_cm^2
  const auto bc   = svec / (-s0);                            // minus the CM velocity
  const auto vc   = toolbox::norm(bc);                       // its magnitude beta_c
  const auto gcm  = sstd::sqrt(q / 2.0);                     // CM Lorentz factor
  const auto vcm  = sstd::sqrt(sstd::max(1.0 - 1.0 / (gcm * gcm), 0.0));
  const auto nm   = toolbox::norm(um);
  const auto omm  = nm > 1e-300 ? um / nm : Vec3(1.0, 0.0, 0.0);  // electron direction
  // A pair at rest as a whole has no boost axis; the electron direction then takes over
  // as the polar axis, and the distribution is symmetric about it (y = 1).
  const auto deg  = vc < 1e-10;
  const auto kvec = deg ? omm : bc / sstd::max(vc, 1e-300);
  auto y = deg
             ? 1.0
             : (gp - gm) / ((gp + gm) * sstd::max(vc, 1e-300) * sstd::max(vcm, 1e-300));
  y      = sstd::clamp(
    y,
    -1.0,
    1.0);  // cosine between the electron's CM direction and the axis
  // azimuth from the (axis, electron) plane
  const auto t = perp_basis(kvec, omm);
  const auto omr =
    cm_direction(y, gcm, vcm, t, rng, ntry);  // photon direction in the CM
  // the two photons, back to back in the CM frame, transformed to the laboratory
  const auto A = boost(bc, gcm, gcm * omr);
  const auto B = boost(bc, gcm, (-gcm) * omr);
  return Products { A.p, B.p };
}

}  // namespace qed::ann
