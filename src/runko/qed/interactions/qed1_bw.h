// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Nonlinear Breit-Wheeler pair production, photon -> electron + positron: the conversion
// rate of a photon in the local field, the closed-form inversion of the pair spectrum, and
// the product sampler. rate() is the rate function r(chi, eps) of a one-body channel on the
// same convention as the emission channel; the threshold XMIN is imposed by the rate itself.
// The field enters through chi_gamma = x b_perp; both leptons leave along the photon, the
// electron taking the fraction r of the quantum parameter (Sarjomaa & Nattila).
// Units: m_e c^2, m_e c.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::bw {

inline constexpr double XMIN   = 2.0;   // photon threshold: one pair's mass
inline constexpr double ZCLIP  = 1e-4;  // deviate clipped off both ends before inversion
inline constexpr double CHI_LO = 0.5;   // validity window of the spectrum shape parameters
inline constexpr double CHI_HI = 1e4;
inline constexpr double TINY   = 1e-300;

/// Auxiliary function T_BW(chi_gamma) of the conversion rate, exponentially small below chi ~ 1.
[[nodiscard]]
constexpr double
  T(double chi_gam)
{
  // (3/16) sqrt(3/2) exp(-8/3chi) in the tunnelling regime, 0.38 chi^(-1/3) once conversion is
  // easy; accurate to 0.6 % from 0.1 to 1e4. Below the floor the exponential has underflowed.
  chi_gam      = sstd::max(chi_gam, 1e-3);
  const auto l = 3.93231893 + sstd::log(chi_gam);
  return 0.37779171 * sstd::exp(-8.0 / (3.0 * chi_gam)) * sstd::pow(chi_gam, -0.333) *
         (1.0 - sstd::exp(-2.003260712099797 * sstd::pow(chi_gam, 0.2688))) *
         (1.0 - 0.6330922 * sstd::exp(-l * l / 29.94554209));
}

/// photon k, field b_perp -> the one-body rate function (chi_gamma/x) T_BW(chi_gamma)
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const double b_perp)
{
  const auto x = toolbox::norm(k);
  if(x < XMIN) { return 0.0; }
  return b_perp * T(x * b_perp);  // chi_gamma/x is b_perp itself
}

struct [[nodiscard]] Shape {
  double a, b, c, n;
};

/// The four Pade parameters of the cumulative probability, at a clamped chi_gamma.
[[nodiscard]]
constexpr Shape
  shape(const double chi_gam)
{
  // Below the floor the approximant acquires a pole (its denominator falls through zero near
  // chi_gamma = 0.3), so the clamp is mandatory; above the ceiling the shape is held.
  const auto x  = sstd::clamp(chi_gam, CHI_LO, CHI_HI);
  const auto y  = x / 100.0;
  const auto n  = 0.25 + 14.3129740 / (1.0 + 250.285446 * (sstd::pow(y, 0.571261959) +
                                                            sstd::pow(y, 2 * 0.571261959))) +
                 0.142285910 * sstd::exp(-4.41463804 * sstd::pow(x / 1000.0, 0.63896860));
  const auto z = x / 150.0;
  const auto a = (0.14549869 * sstd::pow(z, -0.2881441) + 2.10323323 * sstd::pow(z, 0.23634402) - 1.0) *
                 (1.0 + 1.89598777 * sstd::pow(z, -0.74460445));
  const auto b = -1.033388 -
                 (1.0 + 6.74437759 * sstd::pow(y, -1.59943553)) /
                   (1.0 + 2.37682806 * sstd::pow(y, -1.56169871)) -
                 0.13741977 * sstd::pow(x, -2.75389506);
  // c is what makes the approximant pass through p(1/2) = 1/2, so that the half of the
  // spectrum below r = 1/2 can be supplied by reflection
  const auto c = sstd::pow(2.0, -2.0 * n) *
                 (2.0 * a - sstd::pow(2.0, 4.0 * n) - sstd::pow(2.0, 3.0 * n) * b);
  return Shape { a, b, c, n };
}

/// Cumulative probability that the electron takes at most the fraction r of chi_gamma.
[[nodiscard]]
constexpr double
  P(const double chi_gam, const double r)
{
  // Symmetric about r = 1/2: fitted on the upper half, the lower half is reflected.
  const auto up           = r >= 0.5;
  const auto rs           = up ? r : 1.0 - r;
  const auto [a, b, c, n] = shape(chi_gam);
  const auto v            = 1.0 - rs;
  const auto ph = 1.0 - a * sstd::pow(v, 4.0 * n) / (1.0 + b * sstd::pow(v, n) + c * sstd::pow(v, 2.0 * n));
  return up ? ph : 1.0 - ph;
}

/// chi_-/chi_gamma of the produced electron, from a uniform deviate z, by inverting P.
[[nodiscard]]
constexpr double
  r_of(double z, const double chi_gam)
{
  // The deviate is folded into the upper half, the quartic in u = (1 - r)^n solved by Cardano
  // and the depressed quartic, and the result reflected back: the symmetry is exact.
  z                       = sstd::clamp(z, ZCLIP, 1.0 - ZCLIP);
  const auto up           = z >= 0.5;
  const auto zs           = up ? z : 1.0 - z;
  const auto [a, b, c, n] = shape(chi_gam);
  const auto w            = (1.0 - zs) / a;
  const auto wc           = w * c;
  const auto p            = (12.0 * w - wc * wc) / 12.0;
  const auto q   = (2.0 * wc * wc * wc + 72.0 * w * w * c - 27.0 * w * w * b * b) / 216.0;
  const auto dis = sstd::sqrt(sstd::max(q * q / 4.0 + p * p * p / 27.0, 0.0));
  const auto y   = sstd::cbrt(-q / 2.0 + dis) + sstd::cbrt(-q / 2.0 - dis) - wc / 6.0;
  const auto s   = sstd::sqrt(sstd::max(2.0 * y + wc, 0.0));  // the quartic's first square root
  const auto u   = (-s + sstd::sqrt(sstd::max(-2.0 * y + wc - 2.0 * w * b / sstd::max(s, TINY), 0.0))) / 2.0;
  const auto re  = 1.0 - sstd::pow(sstd::max(u, 0.0), 1.0 / n);
  return sstd::clamp(up ? re : 1.0 - re, 0.0, 1.0);
}

struct [[nodiscard]] Products {
  Vec3 um;  // electron
  Vec3 up;  // positron
};

/// photon k, field b_perp -> the electron and positron momenta.
[[nodiscard]]
constexpr Products
  sample(const Vec3& k, const double b_perp, runko::Rng& rng)
{
  // Both leptons leave along the photon, sharing its energy in the ratio the quantum
  // parameters were drawn in, less the pair rest mass: the energies add back to x exactly.
  const auto x  = toolbox::norm(k);
  const auto om = k / sstd::max(x, TINY);
  const auto r  = r_of(rng.uniform(), x * b_perp);
  const auto gm = 1.0 + (x - 2.0) * r;
  const auto gp = 1.0 + (x - 2.0) * (1.0 - r);
  return Products { om * sstd::sqrt(sstd::max(gm * gm - 1.0, 0.0)),
                    om * sstd::sqrt(sstd::max(gp * gp - 1.0, 0.0)) };
}

}  // namespace qed::bw
