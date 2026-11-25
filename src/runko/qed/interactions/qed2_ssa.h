// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed {

inline constexpr double alpha = 7.2973525693e-3;
/// Quantum critical field in G.
inline constexpr double bq_gauss = 4.4140094e13;

/// Returns b_perp = B_perp/B_Q from the magnetic compactness
/// ell_B = sigma_T H U_B/(m_e c^2) = (alpha/3) b_perp^2 (H/lamC).
[[nodiscard]]
constexpr double
  b_from_ellB(const double ell_B, const double size_lamC)
{
  return sstd::sqrt(3.0 * ell_B / (alpha * size_lamC));
}

/// Synchrotron self-absorption, photon + lepton -> lepton.
///
/// The rate is the net of absorption against stimulated emission; it grows without
/// limit towards soft photons, so the bound depends on the field and the photon
/// floor. The field enters only through b_perp = B_perp/B_Q: chi_e = gamma b_perp,
/// chi_gamma = x b_perp. Units: m_e c^2, m_e c, sigma_T.
namespace ssa {

/// Lepton threshold, as for synchrotron emission.
inline constexpr double gamma_min = 3.0;
/// The x/gamma below which the derivative form is used.
inline constexpr double x_rel = 1e-3;

/// Returns the synchrotron function F(nu) = nu int_nu^inf K_5/3.
///
/// F and G(nu) = nu K_2/3(nu) are the fits of Fouka & Ouichaoui, two asymptotes
/// bridged by exp(H1) and 1 - exp(H2); accurate to 0.26 % (F) and 0.035 % (G) over
/// the whole range.
[[nodiscard]]
constexpr double
  F_syn(const double nu)
{
  // F -> f1 nu^(1/3) at nu << 1, F -> f2 e^-nu nu^(1/2) at nu >> 1
  constexpr auto f1 = 2.1495282415344787;  // pi 2^(5/3)/(sqrt 3 Gamma(1/3))
  constexpr auto f2 = 1.2533141373155003;  // sqrt(pi/2)

  const auto x3 = sstd::cbrt(nu);
  const auto H1 = -0.97947838884478688 * nu - 0.83333239129525072 * sstd::sqrt(nu) +
                  0.15541796026816246 * x3;
  const auto H2 = -4.69247165562628882e-2 * nu - 0.70055018056462881 * sstd::sqrt(nu) +
                  1.03876297841949544e-2 * x3;
  return f1 * x3 * sstd::exp(H1) +
         f2 * sstd::exp(-nu) * sstd::sqrt(nu) * (1.0 - sstd::exp(H2));
}

/// Returns K_2/3(nu) = G(nu)/nu.
[[nodiscard]]
constexpr double
  K23(const double nu)
{
  constexpr auto A1c = 0.5 * 1.3541179394264004;  // Gamma(2/3)/2, times (nu/2)^(-2/3)
  constexpr auto A2c = 1.2533141373155003;        // sqrt(pi/2), times nu^(-1/2) e^-nu

  const auto x2 = sstd::sqrt(nu), x3 = sstd::cbrt(nu), x4 = sstd::sqrt(x2);
  const auto H1 = -1.0010216415582440 * nu + 0.88350305221249859 * x2 -
                  3.6240174463901829 * x3 + 0.57393980442916881 * x4;
  const auto H2 = -0.2493940736333195 * nu + 0.9122693061687756 * x2 +
                  1.2051408667145216 * x3 - 5.5227048291651126 * x4;
  return A1c * sstd::pow(nu / 2.0, -2.0 / 3.0) * sstd::exp(H1) +
         A2c / x2 * sstd::exp(-nu) * (1.0 - sstd::exp(H2));
}

/// Returns int_z^inf K_5/3(y) dy.
[[nodiscard]]
constexpr double
  int_K53(const double z)
{
  return F_syn(z) / z;
}

/// Returns the spectral function G(chi_e, chi_gamma) of quantum synchrotron emission;
/// zero for chi_gamma >= chi_e.
[[nodiscard]]
constexpr double
  G_syn(const double chi_e, const double chi_gam)
{
  if(not(chi_gam < chi_e)) { return 0.0; }

  const auto nu =
    sstd::clamp(2.0 * chi_gam / (3.0 * chi_e * (chi_e - chi_gam)), 1e-100, 1e30);
  const auto pref = (sstd::SQRT3 / sstd::TWOPI) * (chi_gam / chi_e);
  return pref * (int_K53(nu) + 1.5 * chi_gam * nu * K23(nu));
}

/// Returns the state-weighted emissivity h(gamma) = u G(gamma b, x b).
[[nodiscard]]
constexpr double
  h_syn(const double g, const double chi_gam, const double b_perp)
{
  return sstd::sqrt(sstd::max(g * g - 1.0, 0.0)) * G_syn(g * b_perp, chi_gam);
}

/// Returns sigma-hat_ssa F/2 for photon k and lepton u, the net of absorption and
/// stimulated emission; zero for gamma < gamma_min or x < x_min.
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const Vec3& u, const double b_perp, const double x_min)
{
  // The opacity is h(gamma + x) - h(gamma); the Moller factor cancels the angular
  // dependence. For x < x_rel gamma the difference is ill-conditioned: x h'(gamma) by a
  // central difference instead, accurate to O(x/gamma).
  const auto x   = toolbox::norm(k);
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  if(gam < gamma_min or x < x_min or x <= 0.0) { return 0.0; }

  const auto chi_gam = x * b_perp;
  const auto us      = sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300));
  double diff;
  if(x < x_rel * gam) {
    const auto d = 1e-6 * gam;
    diff = x * (h_syn(gam + d, chi_gam, b_perp) - h_syn(gam - d, chi_gam, b_perp)) /
           (2.0 * d);
  } else {
    diff = h_syn(gam + x, chi_gam, b_perp) - h_syn(gam, chi_gam, b_perp);
  }

  return (sstd::PI / (8.0 * alpha)) * diff / (us * gam * x * x * x);
}

/// Returns a bound on rate over photons above x_lo and leptons in [g_lo, g_hi].
///
/// It is the maximum over a geometric gamma grid plus the kernel peak
/// gamma^2 = 0.4435 x/b_perp and its neighbours, times 1.05, with x at the floor. The
/// floor maximizes the rate only for gamma b_perp below about 10: above, it rises
/// again as x -> gamma.
[[nodiscard]]
constexpr double
  shmax(const double b_perp, const double x_lo, double g_lo, double g_hi)
{
  g_lo = sstd::max(g_lo, gamma_min);
  g_hi = sstd::max(g_hi, g_lo);

  const auto gpk = sstd::sqrt(0.4435 * x_lo / b_perp);
  const auto k   = Vec3(x_lo, 0.0, 0.0);
  auto r         = 0.0;
  const auto at  = [&](const double g) {
    const auto u = Vec3(sstd::sqrt(g * g - 1.0), 0.0, 0.0);
    r            = sstd::max(r, rate(k, u, b_perp, 0.0));
  };

  constexpr auto n = 16;
  for(auto i = 0; i < n; ++i) {
    at(g_lo * sstd::pow(g_hi / g_lo, static_cast<double>(i) / (n - 1)));
  }
  at(sstd::clamp(gpk / 1.15, g_lo, g_hi));
  at(sstd::clamp(gpk, g_lo, g_hi));
  at(sstd::clamp(gpk * 1.15, g_lo, g_hi));
  return 1.05 * r;
}

/// Returns shmax at the bin edges; 0 for an empty or sub-threshold bin.
[[nodiscard]]
constexpr double
  bound(
    const double b_perp,
    const double x_min,
    const double xlo,
    const double glo,
    const double ghi)
{
  if(not(xlo < 1e300) or not(glo < 1e300) or not(ghi >= gamma_min) or b_perp <= 0.0) {
    return 0.0;
  }
  return shmax(b_perp, sstd::max(xlo, x_min), glo, ghi);
}

/// Returns the absorbing lepton u' with gamma' = gamma + x along u, on the mass shell;
/// the transverse momentum defect is left to the field, as in collinear emission.
[[nodiscard]]
constexpr Vec3
  sample(const Vec3& k, const Vec3& u)
{
  const auto x   = toolbox::norm(k);
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  const auto us  = sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300));
  return u * (sstd::sqrt((gam + x) * (gam + x) - 1.0) / us);
}

/// Returns the largest energy change of the lepton per event, x.
[[nodiscard]]
constexpr double
  dgmax(const double x, const double /*g*/)
{
  return x;
}

}  // namespace ssa
}  // namespace qed
