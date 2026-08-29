// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Quantum synchrotron emission, lepton -> lepton + photon: the photon emission rate in the
// local field, the closed-form inversion of the emission spectrum, and the product sampler.
// rate() is the rate function r(chi, eps) of a one-body channel, the dimensionless factor
// multiplying alpha c/lamC; one-body channels carry no majorant. The field enters through
// b_perp = B_perp/B_Q with chi_e = gamma b_perp; the photon leaves along the lepton momentum
// carrying the fraction r = chi_gamma/chi_e. An energy floor x_min restricts the rate, the
// sampler and the mean to the photons above it, all from the one cumulative probability.
// Everything is elementary (Sarjomaa & Nattila): T_S is a fitted closed form and the
// cumulative of r is a Pade approximant whose inversion is a quartic. Units: m_e c^2, m_e c.

#include "runko/qed/interactions/qed_kinematics.h"

namespace qed::syn {

inline constexpr double GMIN   = 3.0;          // lepton threshold of the channel
inline constexpr double ZCLIP  = 5e-4;         // deviate clipped off both ends before inversion
inline constexpr double CHI_HI = 1e2;          // top of the range the spectrum shape was fitted on
// floor on the quartic's first square root as a divisor
inline constexpr double RMIN   = 1e-9;
// thinnest tail above a floor the sampler can draw from
inline constexpr double FMIN   = 10 * ZCLIP;
inline constexpr double TINY   = 1e-300;

/// Auxiliary function T_S(chi_e) of the emission rate, tending to 5/(2 sqrt 3) as chi_e -> 0.
[[nodiscard]]
constexpr double
  T(double chi_e)
{
  // The photon-number integral of the emission spectrum: the rate function is (chi_e/gamma) T_S.
  // Accurate to 0.5 % from chi_e = 1e-4 to 1e2; the floor keeps the logarithm finite at zero field.
  chi_e         = sstd::max(chi_e, 1e-12);       // zero field floored, not a window
  const auto g1 = 1.0 - 403289.4214570542 * chi_e *
                          sstd::exp(-(chi_e + 0.847336033) * (chi_e + 0.847336033) / 0.0584260036);
  const auto l  = sstd::log(chi_e) + 7.9555351;
  const auto g2 = 1.0 + 3.26342394 * chi_e * sstd::exp(-l * l / 19.18831074);
  return (5.0 / (2.0 * sstd::SQRT3)) *
         (1.0 - sstd::exp(-1.46 * (2.0 * sstd::SQRT3 / 5.0) * sstd::pow(chi_e, -1.0 / 3.0))) *
         g1 * g2;
}

/// Fraction r_min = x_min/x_scale of the quantum parameter below which no photon is emitted.
[[nodiscard]]
constexpr double
  r_floor(const double x_min, const double gam, const bool kinetic)
{
  return x_min / (kinetic ? gam - 1.0 : gam);
}

struct [[nodiscard]] Shape {
  double a, b, n;
};

/// The three Pade parameters of the cumulative probability, at a capped chi_e.
[[nodiscard]]
constexpr Shape
  shape(const double chi_e)
{
  // Downwards the fit extrapolates (the approximant stays a probability to chi_e = 1e-8);
  // upwards the shape is frozen at CHI_HI, the better extrapolation by a wide margin.
  const auto x  = sstd::clamp(chi_e, 1e-12, CHI_HI);  // the floor is zero-field, not a window
  const auto lx = sstd::log(x);
  // the three Pade parameters: amplitude a, denominator coefficient b, exponent n
  const auto a =
    ((1.0 + 0.1648 * sstd::exp(-2.279 * (x - 1.2024) * (x - 1.2024))) * 5.37351138 *
       sstd::pow(x, -0.07620857) / (1.0 + 1.74094714 * sstd::pow(x, 1.29054535)) +
     0.5369956083824988 +
     0.0016994836 * (1.4204061722566785 - 0.5369956083824988) * sstd::pow(x, 1.77253101) /
       (1.0 + 0.00150308200 * sstd::pow(x, 1.77253101))) *
    (1.0 - 0.08494223 * sstd::exp(-4.59296877 * (lx - 1.75725058) * (lx - 1.75725058)));
  const auto b =
    (-1.56564371 +
     (11.3942951 * sstd::pow(x, 0.43384444) + 8.79281269 * sstd::pow(x, 2 * 0.44067827)) *
       sstd::exp(-1.99971294 * sstd::pow(x, 0.46344015)) +
     (-0.38800187 * sstd::pow(x, 2.23683008) + 0.32931219 * sstd::pow(x, 2.35483769)) /
       (1.0 + 0.04777412 * sstd::pow(x, 2.23683008))) *
    sstd::pow(
      0.87786029 + 0.65634725 * sstd::exp(-1.44269312 * (lx + 0.27729158) * (lx + 0.27729158)),
      0.0661953);
  const auto n = 0.25 + 0.926739 * sstd::pow(x, -0.34531647) /
                          (1.0 + 1.26885704 * sstd::pow(x, 0.34531647) +
                           0.30839327 * sstd::pow(x, 1.74193071));
  return Shape { a, b, n };
}

/// Cumulative probability that an emitted photon takes at most the fraction r of chi_e.
[[nodiscard]]
constexpr double
  P(const double chi_e, const double r)
{
  // A Pade form in v = 1 - r^(1/3), anchored at p(0) = 0 by c = a - b - 1 and at p(1) = 1
  // by every term vanishing at r = 1.
  const auto [a, b, n] = shape(chi_e);
  const auto c         = a - b - 1.0;              // anchors p(0) = 0
  const auto v         = 1.0 - sstd::cbrt(r);      // the Pade variable
  return 1.0 - a * sstd::pow(v, 4.0 * n) /
                 (1.0 + b * sstd::pow(v, 2.0 * n) + c * sstd::pow(v, 3.0 * n));
}

/// chi_gamma/chi_e of the emitted photon, from a uniform deviate z, by inverting P.
[[nodiscard]]
constexpr double
  r_of(double z, const double chi_e)
{
  // p = z is a quartic in u = (1 - r^(1/3))^n: the resolvent cubic by Cardano, then one root
  // of the depressed quartic; the sign of hi picks the branch that stays real and positive.
  // The first square root also divides: RMIN floors that divisor alone.
  z                    = sstd::clamp(z, ZCLIP, 1.0 - ZCLIP);
  const auto [a, b, n] = shape(chi_e);
  const auto c         = a - b - 1.0;
  // the quartic a u^4 + B u^2 + C u + D = 0 with z cleared from the denominator ...
  const auto B         = -b * (1.0 - z);
  const auto C         = -c * (1.0 - z);
  const auto D         = z - 1.0;
  // ... depressed to u^4 + K u^2 + M u + N = 0 ...
  const auto K         = B / a - 3.0 * C * C / (8.0 * a * a);
  const auto M         = C * C * C / (8.0 * a * a * a) - B * C / (2.0 * a * a);
  const auto N = D / a + B * C * C / (16.0 * a * a * a) - 3.0 * C * C * C * C / (256.0 * a * a * a * a);
  // ... whose resolvent cubic y^3 + p y + q = 0 Cardano solves
  const auto p   = -(12.0 * N + K * K) / 12.0;
  const auto q   = (-2.0 * K * K * K - 36.0 * K * N + 108.0 * (K * N - M * M / 4.0)) / 216.0;
  const auto dis = sstd::sqrt(sstd::max(q * q / 4.0 + p * p * p / 27.0, 0.0));  // discriminant
  // a real root of the resolvent
  const auto y   = sstd::cbrt(-q / 2.0 + dis) + sstd::cbrt(-q / 2.0 - dis) + K / 6.0;
  // the quartic's roots from the resolvent: the sign of hi picks the real positive branch
  const auto root = sstd::sqrt(sstd::max(2.0 * y - K, 0.0));
  const auto rs   = sstd::max(root, RMIN);        // as a divisor only
  const auto hi   = -2.0 * y - K - 2.0 * M / rs;
  const auto u    = hi >= 0.0
                      ? (root + sstd::sqrt(sstd::max(hi, 0.0))) / 2.0
                      : (-root + sstd::sqrt(sstd::max(-2.0 * y - K + 2.0 * M / rs, 0.0))) / 2.0;
  const auto us   = sstd::max(u, TINY);
  const auto x    = u - C / (4.0 * a);             // undo the depression shift: x = (1 - r^(1/3))^n
  // Classical corner: soft emission with the deviate close to one, where C/(4a) is minute
  // beside u; expanding (u - C/4a)^(1/n) about u keeps the digits that forming x would cancel.
  const auto e    = C / (4.0 * a * n);
  const auto un   = sstd::pow(us, 1.0 / n);
  // second-order expansion of x^(1/n)
  const auto xser = un - un * e / us - e * e * (n - 1.0) * un / (us * us) / 2.0;
  const auto xn = (z < 0.998 or chi_e >= 0.1) ? sstd::pow(sstd::max(x, 0.0), 1.0 / n) : xser;
  const auto r  = x < 0.0 ? 1.0 : (1.0 - xn) * (1.0 - xn) * (1.0 - xn);  // negative x: round-off
  return sstd::clamp(r, 0.0, 1.0);
}

/// Fraction of the emitted photons above the floor r_min: 1 - P(chi_e, r_min), or 0 when the
/// tail is thinner than the clip the inversion applies (the channel is closed for that lepton).
[[nodiscard]]
constexpr double
  frac(const double chi_e, const double r_min)
{
  const auto f = 1.0 - P(chi_e, r_min);
  return f >= FMIN ? f : 0.0;
}

/// lepton u, field b_perp -> the one-body rate function (chi_e/gamma) T_S(chi_e), counting
/// only the photons above the floor x_min.
[[nodiscard]]
constexpr double
  rate(const Vec3& u, const double b_perp, const double x_min, const bool kinetic)
{
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  if(gam < GMIN) { return 0.0; }
  const auto chi_e = gam * b_perp;
  return (chi_e / gam) * T(chi_e) * frac(chi_e, r_floor(x_min, gam, kinetic));
}

/// Quantum correction of the radiated power, falling from g = 1 in the classical limit.
[[nodiscard]]
constexpr double
  g(const double chi_e)
{
  return sstd::pow(1.0 + 4.8 * (1.0 + chi_e) * sstd::log(1.0 + 1.7 * chi_e) + 2.44 * chi_e * chi_e,
                   -2.0 / 3.0);
}

/// Mean fraction <chi_gamma>/chi_e that one drawn photon carries, among the photons above r_min.
[[nodiscard]]
constexpr double
  rmean(double chi_e, const double r_min)
{
  // The mean of what the sampler draws, from the tail integral of its own cumulative:
  // int_{r_min}^1 r dP = r_min (1 - P0) + int_{r_min}^1 (1 - P) dr, Simpson in ln r from a
  // few decades below the photon energy scale r ~ chi_e; the sliver below the first node is
  // added in closed form.
  chi_e           = sstd::max(chi_e, 1e-12);       // zero field floored, as in shape()
  const auto P0   = P(chi_e, r_min);               // cumulative at the floor
  // first node, below the spectrum's scale
  auto r_lo       = sstd::max(r_min, 1e-7 * sstd::min(chi_e, 1.0));
  r_lo            = sstd::min(r_lo, 1.0 - 1e-12);
  constexpr int n = 96;                            // Simpson intervals in ln r
  const auto lr0  = sstd::log(r_lo);
  auto sum        = 0.0;
  // Simpson sum of (1 - P) r over ln r from r_lo to 1
  for(auto i = 0; i <= n; ++i) {
    const auto r = sstd::exp(lr0 * (1.0 - static_cast<double>(i) / n));
    const auto w = (i == 0 or i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    sum += w * (1.0 - P(chi_e, r)) * r;
  }
  const auto h    = -lr0 / n;
  const auto tail = h / 3.0 * sum + (r_lo - r_min) * (1.0 - P(chi_e, r_lo));
  return r_min + tail / sstd::max(1.0 - P0, TINY);
}

/// Mean energy of a drawn photon above the floor: the scale of the multiplicity of a merged emission.
[[nodiscard]]
constexpr double
  xmean(const double chi_e, const double gam, const bool kinetic, const double x_min)
{
  const auto scale = kinetic ? gam - 1.0 : gam;
  return scale * rmean(chi_e, r_floor(x_min, gam, kinetic));
}

struct [[nodiscard]] Products {
  Vec3 k;   // emitted photon
  Vec3 u1;  // recoiled lepton
};

/// lepton u, field b_perp -> the emitted photon (above the floor x_min) and the recoiled lepton.
[[nodiscard]]
constexpr Products
  sample(const Vec3& u, const double b_perp, runko::Rng& rng, const bool kinetic, const double x_min)
{
  // Collinear emission: the photon leaves along the lepton momentum and the lepton keeps its
  // direction while its energy drops by the photon energy, rescaled on the mass shell -- the
  // exact inverse of the self-absorption product. The floor restricts the deviate to the part
  // of the cumulative above it, at no rejection.
  const auto gam   = sstd::sqrt(1.0 + toolbox::dot(u, u));  // Lorentz factor of the emitter
  const auto us    = sstd::sqrt(sstd::max(gam * gam - 1.0, TINY));  // its momentum, to be rescaled
  const auto om    = u / us;                       // emission direction, along the momentum
  const auto chi_e = gam * b_perp;
  // cumulative at the floor, 0 without one
  const auto P0    = P(chi_e, r_floor(x_min, gam, kinetic));
  // fraction drawn above the floor
  const auto r     = r_of(P0 + (1.0 - P0) * rng.uniform(), chi_e);
  const auto x     = (kinetic ? gam - 1.0 : gam) * r;  // photon energy
  const auto gam1  = sstd::max(gam - x, 1.0);      // lepton energy after the emission
  return Products { x * om, om * sstd::sqrt(sstd::max(gam1 * gam1 - 1.0, 0.0)) };
}

}  // namespace qed::syn
