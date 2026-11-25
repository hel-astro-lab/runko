// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Nonlinear Breit-Wheeler pair production, photon -> electron + positron, in the local
// field: the one-body rate function r(chi, eps) (zero below x_threshold), the closed-form
// inversion of the pair spectrum (Sarjomaa & Nattila), and the product sampler. The
// field enters through chi_gamma = x b_perp. Units: m_e c^2, m_e c.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::bw {

inline constexpr double x_threshold  = 2.0;   // photon threshold: one pair's mass
inline constexpr double z_clip = 1e-4;  // deviate clipped off both ends before inversion
inline constexpr double chi_lo = 0.5;  // validity window of the spectrum shape
inline constexpr double chi_hi = 1e4;
inline constexpr double tiny   = 1e-300;

/// Returns T_BW(chi_gamma) of the conversion rate: (3/16) sqrt(3/2) exp(-8/(3 chi)) in
/// the tunnelling regime, 0.38 chi^(-1/3) at large chi; fit accurate to 0.6 % from 0.1
/// to 1e4.
[[nodiscard]]
constexpr double
  T(double chi_gam)
{
  chi_gam      = sstd::max(chi_gam, 1e-3);  // below, the exponential underflows anyway
  const auto l = 3.93231893 + sstd::log(chi_gam);
  return 0.37779171 * sstd::exp(-8.0 / (3.0 * chi_gam)) * sstd::pow(chi_gam, -0.333) *
         (1.0 - sstd::exp(-2.003260712099797 * sstd::pow(chi_gam, 0.2688))) *
         (1.0 - 0.6330922 * sstd::exp(-l * l / 29.94554209));
}

/// Returns the one-body rate function (chi_gamma/x) T_BW(chi_gamma) = b_perp T_BW.
[[nodiscard]]
constexpr double
  rate(const Vec3& k, const double b_perp)
{
  const auto x = toolbox::norm(k);
  if(x < x_threshold) { return 0.0; }
  return b_perp * T(x * b_perp);
}

/// Pade parameters of the cumulative pair spectrum: amplitude a, denominator
/// coefficients b and c, exponent n.
struct Shape {
  double a, b, c, n;
};

/// Returns the Pade parameters at chi_gamma clamped to [chi_lo, chi_hi]. The
/// denominator crosses zero at chi = 0.205, so the lower clamp is mandatory.
[[nodiscard]]
constexpr Shape
  shape(const double chi_gam)
{
  const auto x = sstd::clamp(chi_gam, chi_lo, chi_hi);
  const auto y = x / 100.0;

  const auto n =
    0.25 +
    14.3129740 /
      (1.0 + 250.285446 * (sstd::pow(y, 0.571261959) + sstd::pow(y, 2 * 0.571261959))) +
    0.142285910 * sstd::exp(-4.41463804 * sstd::pow(x / 1000.0, 0.63896860));
  const auto z = x / 150.0;
  const auto a = (0.14549869 * sstd::pow(z, -0.2881441) +
                  2.10323323 * sstd::pow(z, 0.23634402) - 1.0) *
                 (1.0 + 1.89598777 * sstd::pow(z, -0.74460445));
  const auto b = -1.033388 -
                 (1.0 + 6.74437759 * sstd::pow(y, -1.59943553)) /
                   (1.0 + 2.37682806 * sstd::pow(y, -1.56169871)) -
                 0.13741977 * sstd::pow(x, -2.75389506);

  // c makes p(1/2) = 1/2, so the half r < 1/2 follows by reflection
  const auto c = sstd::pow(2.0, -2.0 * n) *
                 (2.0 * a - sstd::pow(2.0, 4.0 * n) - sstd::pow(2.0, 3.0 * n) * b);
  return Shape { a, b, c, n };
}

/// Returns the probability that the electron takes at most the fraction r of chi_gamma:
/// 1 - a v^(4n)/(1 + b v^n + c v^(2n)), v = 1 - r, on r >= 1/2; reflected below.
[[nodiscard]]
constexpr double
  P(const double chi_gam, const double r)
{
  const auto up           = r >= 0.5;
  const auto rs           = up ? r : 1.0 - r;
  const auto [a, b, c, n] = shape(chi_gam);
  const auto v            = 1.0 - rs;
  const auto ph = 1.0 - a * sstd::pow(v, 4.0 * n) /
                          (1.0 + b * sstd::pow(v, n) + c * sstd::pow(v, 2.0 * n));
  return up ? ph : 1.0 - ph;
}

/// Returns chi_-/chi_gamma of the electron from a uniform deviate z by inverting P.
[[nodiscard]]
constexpr double
  r_of(double z, const double chi_gam)
{
  // fold z into the upper half, invert there, reflect back
  z                       = sstd::clamp(z, z_clip, 1.0 - z_clip);
  const auto up           = z >= 0.5;
  const auto zs           = up ? z : 1.0 - z;
  const auto [a, b, c, n] = shape(chi_gam);

  // quartic u^4 - w c u^2 - w b u - w = 0 in u = (1 - r)^n, w = (1 - z)/a; its
  // resolvent cubic depressed to t^3 + p t + q = 0, y = t - wc/6, solved by Cardano
  const auto w  = (1.0 - zs) / a;
  const auto wc = w * c;
  const auto p  = (12.0 * w - wc * wc) / 12.0;
  const auto q = (2.0 * wc * wc * wc + 72.0 * w * w * c - 27.0 * w * w * b * b) / 216.0;
  const auto dis = sstd::sqrt(sstd::max(q * q / 4.0 + p * p * p / 27.0, 0.0));
  const auto y   = sstd::cbrt(-q / 2.0 + dis) + sstd::cbrt(-q / 2.0 - dis) - wc / 6.0;

  // the real positive root of the quartic
  const auto s = sstd::sqrt(sstd::max(2.0 * y + wc, 0.0));
  const auto u =
    (-s +
     sstd::sqrt(sstd::max(-2.0 * y + wc - 2.0 * w * b / sstd::max(s, tiny), 0.0))) /
    2.0;

  const auto re = 1.0 - sstd::pow(sstd::max(u, 0.0), 1.0 / n);
  return sstd::clamp(up ? re : 1.0 - re, 0.0, 1.0);
}

struct Products {
  Vec3 um;  // electron
  Vec3 up;  // positron
};

/// Returns the electron and positron of photon k in the field b_perp: both along the
/// photon, sharing x - 2 in the drawn ratio r, so the energies add to x exactly.
[[nodiscard]]
constexpr Products
  sample(const Vec3& k, const double b_perp, runko::Rng& rng)
{
  const auto x  = toolbox::norm(k);
  const auto om = k / sstd::max(x, tiny);
  const auto r  = r_of(rng.uniform(), x * b_perp);

  const auto gm = 1.0 + (x - 2.0) * r;
  const auto gp = 1.0 + (x - 2.0) * (1.0 - r);
  return Products { om * sstd::sqrt(sstd::max(gm * gm - 1.0, 0.0)),
                    om * sstd::sqrt(sstd::max(gp * gp - 1.0, 0.0)) };
}

}  // namespace qed::bw
