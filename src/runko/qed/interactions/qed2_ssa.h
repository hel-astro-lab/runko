// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Synchrotron self-absorption, photon + lepton -> lepton: the net rate of absorption against
// stimulated emission, its field-dependent bound, and the single product the conservation
// laws fix. rate() is sigma-hat_ssa F/2 on the convention of the other two-body channels,
// but it admits no constant bound: it grows without limit towards soft photons, so shmax()
// is a function of the field and of the low-energy floor. The field enters only through
// b_perp = B_perp/B_Q: chi_e = gamma b_perp, chi_gamma = x b_perp.
// Units: m_e c^2, m_e c, sigma_T.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed {

inline constexpr double ALPHA = 7.2973525693e-3;  // fine structure constant
inline constexpr double BQ_G  = 4.4140094e13;     // quantum critical field, G

/// b_perp = B_perp/B_Q from the magnetic compactness ell_B = sigma_T H U_B/(m_e c^2):
/// ell_B = (alpha/3) b_perp^2 (H/lamC), so the size of the system is an independent scale.
[[nodiscard]]
constexpr double
  b_from_ellB(const double ell_B, const double size_lamC)
{
  return sstd::sqrt(3.0 * ell_B / (ALPHA * size_lamC));
}

namespace ssa {

inline constexpr double GMIN = 3.0;   // lepton threshold, as for synchrotron emission
inline constexpr double XREL = 1e-3;  // x/gamma below which the rate takes its derivative form

// ---------------------------------------------------------------- synchrotron functions
// F(nu) = nu int_nu^inf K_5/3 and G(nu) = nu K_2/3(nu) as the elementary fits of Fouka &
// Ouichaoui: each is its two asymptotes bridged by exp(H1) and 1 - exp(H2), H a short sum of
// fractional powers, accurate to 0.26 % (F) and 0.035 % (G) over the whole range.

[[nodiscard]]
constexpr double
  F_syn(const double nu)
{
  // pi 2^(5/3)/(sqrt 3 Gamma(1/3)): F -> F1 nu^(1/3) at nu << 1
  constexpr auto F1 = 2.1495282415344787;
  constexpr auto F2 = 1.2533141373155003;  // sqrt(pi/2): F -> F2 e^-nu nu^(1/2) at nu >> 1
  const auto x3     = sstd::cbrt(nu);      // nu^(1/3)
  const auto H1     = -0.97947838884478688 * nu - 0.83333239129525072 * sstd::sqrt(nu) +
                      0.15541796026816246 * x3;  // exp(H1) bridges away from the soft asymptote
  const auto H2 = -4.69247165562628882e-2 * nu - 0.70055018056462881 * sstd::sqrt(nu) +
                  1.03876297841949544e-2 * x3;  // 1 - exp(H2) bridges into the hard one
  return F1 * x3 * sstd::exp(H1) + F2 * sstd::exp(-nu) * sstd::sqrt(nu) * (1.0 - sstd::exp(H2));
}

[[nodiscard]]
constexpr double
  K23(const double nu)
{
  // K_2/3 itself; G(nu) = nu K_2/3(nu)
  constexpr auto A1c = 0.5 * 1.3541179394264004;  // Gamma(2/3)/2, times (nu/2)^(-2/3)
  constexpr auto A2c = 1.2533141373155003;        // sqrt(pi/2), times nu^(-1/2) e^-nu
  // nu^(1/2), nu^(1/3), nu^(1/4)
  const auto x2 = sstd::sqrt(nu), x3 = sstd::cbrt(nu), x4 = sstd::sqrt(x2);
  const auto H1 = -1.0010216415582440 * nu + 0.88350305221249859 * x2 -
                  3.6240174463901829 * x3 + 0.57393980442916881 * x4;
  const auto H2 = -0.2493940736333195 * nu + 0.9122693061687756 * x2 +
                  1.2051408667145216 * x3 - 5.5227048291651126 * x4;
  return A1c * sstd::pow(nu / 2.0, -2.0 / 3.0) * sstd::exp(H1) +
         A2c / x2 * sstd::exp(-nu) * (1.0 - sstd::exp(H2));
}

/// int_z^inf K_5/3(y) dy
[[nodiscard]]
constexpr double
  int_K53(const double z)
{
  return F_syn(z) / z;
}

/// Spectral function G(chi_e, chi_gamma) of the quantum synchrotron spectrum.
[[nodiscard]]
constexpr double
  G_syn(const double chi_e, const double chi_gam)
{
  if(not(chi_gam < chi_e)) { return 0.0; }  // a photon cannot carry more than the lepton has
  // argument of the synchrotron functions
  const auto nu   = sstd::clamp(2.0 * chi_gam / (3.0 * chi_e * (chi_e - chi_gam)), 1e-100, 1e30);
  // sqrt(3)/(2 pi) chi_gamma/chi_e
  const auto pref = (1.7320508075688772 / 6.283185307179586) * (chi_gam / chi_e);
  return pref * (int_K53(nu) + 1.5 * chi_gam * nu * K23(nu));
}

/// state-weighted emissivity h(gamma) = u G(gamma b, x b) of the absorption step
[[nodiscard]]
constexpr double
  h_syn(const double g, const double chi_gam, const double b_perp)
{
  return sstd::sqrt(sstd::max(g * g - 1.0, 0.0)) * G_syn(g * b_perp, chi_gam);
}

/// photon k, lepton u -> sigma-hat_ssa F/2, the net of absorption and stimulated emission.
/// Zero below the thresholds (gamma < GMIN, x < x_min).
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const Vec3& u, const double b_perp, const double x_min)
{
  // The two differ only at order x/gamma and the opacity is their difference, so the channel
  // carries the difference of h(gamma) across the absorption step gamma -> gamma + x. The
  // Moller factor cancels the angular dependence, so the rate depends on x and gamma alone.
  // That difference is ill-conditioned for a soft photon: below x = XREL gamma the derivative
  // h'(gamma), a central difference of the same smooth fit, is accurate to (x/gamma)^2 there.
  const auto x   = toolbox::norm(k);
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  if(gam < GMIN or x < x_min or x <= 0.0) { return 0.0; }
  const auto chi_gam = x * b_perp;                         // the same for both states of the lepton
  const auto us      = sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300));  // lepton momentum
  double diff;                                             // h(gamma + x) - h(gamma)
  if(x < XREL * gam) {
    const auto d = 1e-6 * gam;
    diff = x * (h_syn(gam + d, chi_gam, b_perp) - h_syn(gam - d, chi_gam, b_perp)) / (2.0 * d);
  } else {
    diff = h_syn(gam + x, chi_gam, b_perp) - h_syn(gam, chi_gam, b_perp);
  }
  return (3.141592653589793 / (8.0 * ALPHA)) * diff / (us * gam * x * x * x);
}

/// Bound on rate over photons above x_lo and leptons in [g_lo, g_hi], by direct maximization:
/// the rate falls with x, so the bound sits at the floor and only gamma is scanned, on a
/// geometric grid with the peak gamma^2 = 0.4435 x/b_perp of the kernel and its neighbours.
[[nodiscard]]
constexpr double
  shmax(const double b_perp, const double x_lo, double g_lo, double g_hi)
{
  g_lo           = sstd::max(g_lo, GMIN);
  g_hi           = sstd::max(g_hi, g_lo);
  const auto gpk = sstd::sqrt(0.4435 * x_lo / b_perp);    // gamma of the kernel's peak at the floor
  const auto k   = Vec3(x_lo, 0.0, 0.0);                  // a photon at the floor
  auto r         = 0.0;                                   // running maximum of the rate
  const auto at  = [&](const double g) {                  // the rate at gamma, head-on
    const auto u = Vec3(sstd::sqrt(g * g - 1.0), 0.0, 0.0);
    r            = sstd::max(r, rate(k, u, b_perp, 0.0));
  };
  constexpr auto n = 16;
  // geometric scan of the gamma range, then the peak and its neighbours
  for(auto i = 0; i < n; ++i) {
    at(g_lo * sstd::pow(g_hi / g_lo, static_cast<double>(i) / (n - 1)));
  }
  at(sstd::clamp(gpk / 1.15, g_lo, g_hi));
  at(sstd::clamp(gpk, g_lo, g_hi));
  at(sstd::clamp(gpk * 1.15, g_lo, g_hi));
  return 1.05 * r;  // margin over the grid maximum
}

/// Bound over a bin pair from the bin edges; 0 for an empty or sub-threshold bin.
[[nodiscard]]
constexpr double
  bound(const double b_perp, const double x_min, const double xlo, const double glo, const double ghi)
{
  if(not(xlo < 1e300) or not(glo < 1e300) or not(ghi >= GMIN)) { return 0.0; }
  return shmax(b_perp, sstd::max(xlo, x_min), glo, ghi);
}

/// photon k, lepton u -> the single product u', with gamma' = gamma + x: the lepton keeps its
/// direction and takes the photon energy, on the mass shell; the transverse momentum defect
/// is left to the field, as in the collinear emission this is the inverse of.
[[nodiscard]]
constexpr Vec3
  sample(const Vec3& k, const Vec3& u)
{
  const auto x   = toolbox::norm(k);                      // photon energy, all absorbed
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));  // Lorentz factor of the absorbing lepton
  const auto us  = sstd::sqrt(sstd::max(gam * gam - 1.0, 1e-300));  // its momentum, to be rescaled
  return u * (sstd::sqrt((gam + x) * (gam + x) - 1.0) / us);
}

/// Largest energy change of the lepton per event, for the accumulation factor.
[[nodiscard]]
constexpr double
  dgmax(const double x, const double /*g*/)
{
  return x;
}

}  // namespace ssa
}  // namespace qed
