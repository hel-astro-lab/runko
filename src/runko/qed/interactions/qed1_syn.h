// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/qed/interactions/qed_kinematics.h"

#include <array>

/// Quantum synchrotron emission, lepton -> lepton + photon, in the local field.
///
/// The one-body rate function (the factor of alpha c/lamC), the inversion of the
/// emission spectrum, the mean photon energy above a floor x_min and the product
/// sampler are all in closed form (Sarjomaa & Nattila). The quantum parameter is
/// chi_e = gamma b_perp with b_perp = B_perp/B_Q; the photon takes the fraction
/// r = chi_gamma/chi_e. Units: m_e c^2, m_e c.
namespace qed::syn {

/// Lepton threshold of the channel.
inline constexpr double gamma_min = 3.0;
/// Deviate clipped off both ends.
inline constexpr double z_clip = 5e-4;
/// Top of the fitted range of the shape.
inline constexpr double chi_hi = 1e2;
/// Floor of a divisor in the inversion.
inline constexpr double r_tiny = 1e-9;
/// Thinnest drawable tail.
inline constexpr double frac_min = 10 * z_clip;
inline constexpr double tiny     = 1e-300;

/// Returns T_S(chi_e), the photon-number integral of the spectrum.
///
/// The rate function is (chi_e/gamma) T_S. T_S tends to 5/(2 sqrt 3) as chi_e -> 0;
/// the fit is accurate to 0.5 % from chi_e = 1e-4 to 1e2.
[[nodiscard]]
constexpr double
  T(double chi_e)
{
  chi_e = sstd::max(chi_e, 1e-12);  // zero field: keeps the logarithm finite

  const auto g1 =
    1.0 - 403289.4214570542 * chi_e *
            sstd::exp(-(chi_e + 0.847336033) * (chi_e + 0.847336033) / 0.0584260036);
  const auto l  = sstd::log(chi_e) + 7.9555351;
  const auto g2 = 1.0 + 3.26342394 * chi_e * sstd::exp(-l * l / 19.18831074);
  return (5.0 / (2.0 * sstd::SQRT3)) *
         (1.0 -
          sstd::exp(-1.46 * (2.0 * sstd::SQRT3 / 5.0) * sstd::pow(chi_e, -1.0 / 3.0))) *
         g1 * g2;
}

/// Returns r_min = x_min/x_scale, x_scale = gamma - 1 (kinetic) or gamma.
[[nodiscard]]
constexpr double
  r_floor(const double x_min, const double gam, const bool kinetic)
{ return x_min / (kinetic ? gam - 1.0 : gam); }

/// Pade parameters (a, b, n) of the cumulative spectrum and the parameters
/// (l0, linf, kappa) of the effective point of its mean (see rmean).
struct Shape {
  double a, b, n, l0, linf, kappa;
};

// Effective point of the mean: l0, linf, ln kappa as a cubic in
// s = (ln chi_e - s_mid)/s_half plus one logistic step, fitted on 1e-12 <= chi_e <= 10.
inline constexpr double s_mid = -12.66422, s_half = 14.9668;
inline constexpr std::array<double, 7> l0_fit = { 4.749836e-01, 5.957273e-02, 2.173267e-02,
                                              1.045624e-02, 2.865377e-01, 5.673384e-01,
                                              1.583990e+00 };
inline constexpr std::array<double, 7> linf_fit    = { 8.268931e-01,  -1.466373e-02,
                                                   -1.051902e-02, -8.028332e-03,
                                                   8.251289e-02,  6.572695e-02,
                                                   7.165592e-01 };
inline constexpr std::array<double, 7> lnkappa_fit = { -9.414223e-01, 3.268251e-01,
                                                   6.384624e-02,  2.694743e-02,
                                                   2.905601e+00,  2.529502e+00,
                                                   2.155030e+00 };
/// Stretching exponent of the point.
inline constexpr double p_stretch = 1.060691e+00;

/// Returns c0 + c1 s + c2 s^2 + c3 s^3 + E/(1 + exp(-(t - t_s)/w_s)) at t = ln chi_e.
[[nodiscard]]
constexpr double
  curve(const std::array<double, 7>& q, const double t)
{
  const auto s = (t - s_mid) / s_half;
  return q[0] + q[1] * s + q[2] * s * s + q[3] * s * s * s +
         q[4] / (1.0 + sstd::exp(-(t - q[5]) / q[6]));
}

/// Returns the shape at chi_e.
///
/// Downwards the fit extrapolates (still a probability to chi_e = 1e-8; 1e-12 floors
/// zero field); above chi_hi it is frozen, the far better extrapolation.
[[nodiscard]]
constexpr Shape
  shape(const double chi_e)
{
  const auto x  = sstd::clamp(chi_e, 1e-12, chi_hi);
  const auto lx = sstd::log(x);

  const auto a =
    ((1.0 + 0.1648 * sstd::exp(-2.279 * (x - 1.2024) * (x - 1.2024))) * 5.37351138 *
       sstd::pow(x, -0.07620857) / (1.0 + 1.74094714 * sstd::pow(x, 1.29054535)) +
     0.5369956083824988 +
     0.0016994836 * (1.4204061722566785 - 0.5369956083824988) *
       sstd::pow(x, 1.77253101) / (1.0 + 0.00150308200 * sstd::pow(x, 1.77253101))) *
    (1.0 - 0.08494223 * sstd::exp(-4.59296877 * (lx - 1.75725058) * (lx - 1.75725058)));
  const auto b =
    (-1.56564371 +
     (11.3942951 * sstd::pow(x, 0.43384444) +
      8.79281269 * sstd::pow(x, 2 * 0.44067827)) *
       sstd::exp(-1.99971294 * sstd::pow(x, 0.46344015)) +
     (-0.38800187 * sstd::pow(x, 2.23683008) + 0.32931219 * sstd::pow(x, 2.35483769)) /
       (1.0 + 0.04777412 * sstd::pow(x, 2.23683008))) *
    sstd::pow(
      0.87786029 +
        0.65634725 * sstd::exp(-1.44269312 * (lx + 0.27729158) * (lx + 0.27729158)),
      0.0661953);
  const auto n = 0.25 + 0.926739 * sstd::pow(x, -0.34531647) /
                          (1.0 + 1.26885704 * sstd::pow(x, 0.34531647) +
                           0.30839327 * sstd::pow(x, 1.74193071));

  const auto t = sstd::log(x);
  return Shape { a, b, n, curve(l0_fit, t), curve(linf_fit, t), sstd::exp(curve(lnkappa_fit, t)) };
}

/// Returns the probability that a photon takes at most the fraction r of chi_e.
///
/// It is 1 - a u^4/(1 + b u^2 + c u^3), u = (1 - r^(1/3))^n, with c = a - b - 1 so
/// that p(0) = 0.
[[nodiscard]]
constexpr double
  P(const Shape& s, const double r)
{
  const auto c  = s.a - s.b - 1.0;
  const auto v  = 1.0 - sstd::cbrt(r);
  const auto vn = sstd::pow(v, s.n);
  return 1.0 - s.a * vn * vn * vn * vn / (1.0 + s.b * vn * vn + c * vn * vn * vn);
}

[[nodiscard]]
constexpr double
  P(const double chi_e, const double r)
{ return P(shape(chi_e), r); }

/// Returns chi_gamma/chi_e of the photon from a uniform deviate z by inverting P.
[[nodiscard]]
constexpr double
  r_of(double z, const Shape& s, const double chi_e)
{
  z            = sstd::clamp(z, z_clip, 1.0 - z_clip);
  const auto a = s.a, b = s.b, n = s.n;
  const auto c = a - b - 1.0;

  // P = z as the quartic a u^4 + B u^2 + C u + D = 0, depressed to
  // u^4 + K u^2 + M u + N = 0, whose resolvent cubic y^3 + p y + q = 0 Cardano solves
  const auto B = -b * (1.0 - z);
  const auto C = -c * (1.0 - z);
  const auto D = z - 1.0;
  const auto K = B / a - 3.0 * C * C / (8.0 * a * a);
  const auto M = C * C * C / (8.0 * a * a * a) - B * C / (2.0 * a * a);
  const auto N = D / a + B * C * C / (16.0 * a * a * a) -
                 3.0 * C * C * C * C / (256.0 * a * a * a * a);
  const auto p = -(12.0 * N + K * K) / 12.0;
  const auto q =
    (-2.0 * K * K * K - 36.0 * K * N + 108.0 * (K * N - M * M / 4.0)) / 216.0;
  const auto dis = sstd::sqrt(sstd::max(q * q / 4.0 + p * p * p / 27.0, 0.0));
  const auto y   = sstd::cbrt(-q / 2.0 + dis) + sstd::cbrt(-q / 2.0 - dis) + K / 6.0;

  // the sign of hi picks the branch that stays real and positive
  const auto root = sstd::sqrt(sstd::max(2.0 * y - K, 0.0));
  const auto rs   = sstd::max(root, r_tiny);
  const auto hi   = -2.0 * y - K - 2.0 * M / rs;
  const auto u =
    hi >= 0.0 ? (root + sstd::sqrt(sstd::max(hi, 0.0))) / 2.0
              : (-root + sstd::sqrt(sstd::max(-2.0 * y - K + 2.0 * M / rs, 0.0))) / 2.0;
  const auto us = sstd::max(u, tiny);
  const auto x  = u - C / (4.0 * a);  // (1 - r^(1/3))^n, the depression shift undone

  // Classical corner (soft emission, z near 1): C/(4a) is minute beside u, so x^(1/n)
  // is expanded about u to second order instead of formed by cancellation
  const auto e    = C / (4.0 * a * n);
  const auto un   = sstd::pow(us, 1.0 / n);
  const auto xser = un - un * e / us - e * e * (n - 1.0) * un / (us * us) / 2.0;
  const auto xn =
    (z < 0.998 or chi_e >= 0.1) ? sstd::pow(sstd::max(x, 0.0), 1.0 / n) : xser;

  // x < 0 is round-off
  const auto r = x < 0.0 ? 1.0 : (1.0 - xn) * (1.0 - xn) * (1.0 - xn);
  return sstd::clamp(r, 0.0, 1.0);
}

[[nodiscard]]
constexpr double
  r_of(const double z, const double chi_e)
{ return r_of(z, shape(chi_e), chi_e); }

/// Returns the fraction 1 - P(chi_e, r_min) of photons above the floor, or 0 (channel
/// closed) when it is thinner than frac_min.
[[nodiscard]]
constexpr double
  frac(const Shape& s, const double r_min)
{
  const auto f = 1.0 - P(s, r_min);
  return f >= frac_min ? f : 0.0;
}

[[nodiscard]]
constexpr double
  frac(const double chi_e, const double r_min)
{ return frac(shape(chi_e), r_min); }

/// Returns the one-body rate function (chi_e/gamma) T_S(chi_e) of the photons above
/// x_min, for a lepton of Lorentz factor gam at chi_e with shape s.
[[nodiscard]]
constexpr double
  rate(
    const Shape& s,
    const double chi_e,
    const double gam,
    const double x_min,
    const bool kinetic)
{
  if(gam < gamma_min) { return 0.0; }
  return (chi_e / gam) * T(chi_e) * frac(s, r_floor(x_min, gam, kinetic));
}

/// Returns the same for lepton u in the field b_perp.
[[nodiscard]]
constexpr double
  rate(const Vec3& u, const double b_perp, const double x_min, const bool kinetic)
{
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  if(gam < gamma_min) { return 0.0; }

  const auto chi_e = gam * b_perp;
  return rate(shape(chi_e), chi_e, gam, x_min, kinetic);
}

/// Returns the quantum correction g(chi_e) of the radiated power; g = 1 classically.
[[nodiscard]]
constexpr double
  g(const double chi_e)
{
  return sstd::pow(
    1.0 + 4.8 * (1.0 + chi_e) * sstd::log(1.0 + 1.7 * chi_e) + 2.44 * chi_e * chi_e,
    -2.0 / 3.0);
}

/// Returns the mean <chi_gamma>/chi_e of the photons the sampler draws above r_min, in
/// closed form.
///
/// Multiplicities of merged emissions divide by this drawn mean, not by the exact
/// spectrum's.
[[nodiscard]]
constexpr double
  rmean(const Shape& s, const double r_min)
{
  // <r> = r_min + int_{r_min}^1 (1 - P) dr / (1 - P0). For the numerator a u^4 alone,
  // with rho = r_min^(1/3) and m = 4n, the integral is exactly
  //   Delta_1 = 3 (1 - rho) [2 + 2 rho (m+1) + rho^2 (m+1)(m+2)] / ((m+1)(m+2)(m+3)).
  // The denominator D(u) = 1 + b u^2 + c u^3 contributes D(u_m)/D(lambda u_m); the
  // effective point lambda relaxes from l0 (no floor) to linf deep in the tail as
  // exp(-kappa z^p), z = -4 ln u_m the tail exponent. Within 0.3 % of the moment
  // for chi_e <= 1e-2, 0.5 % to 1, 0.8 % to 10 and 2 % extrapolated to 100, over
  // every floor that leaves the channel open.
  const auto c   = s.a - s.b - 1.0;
  const auto m   = 4.0 * s.n;
  const auto rho = sstd::cbrt(r_min);
  const auto v   = sstd::max(1.0 - rho, tiny);
  const auto d1  = 3.0 * v *
                   (2.0 + 2.0 * rho * (m + 1.0) + rho * rho * (m + 1.0) * (m + 2.0)) /
                   ((m + 1.0) * (m + 2.0) * (m + 3.0));

  const auto z = -m * sstd::log(v);
  const auto u = sstd::exp(-z / 4.0);  // u_m = v_m^n
  const auto lam =
    s.linf - (s.linf - s.l0) * sstd::exp(-s.kappa * sstd::pow(z, p_stretch));
  const auto D = [&](const double t) { return 1.0 + s.b * t * t + c * t * t * t; };
  return r_min + d1 * D(u) / D(lam * u);
}

[[nodiscard]]
constexpr double
  rmean(const double chi_e, const double r_min)
{ return rmean(shape(chi_e), r_min); }

/// Returns the mean energy of a drawn photon above x_min, x_scale rmean.
[[nodiscard]]
constexpr double
  xmean(const Shape& s, const double gam, const bool kinetic, const double x_min)
{
  const auto scale = kinetic ? gam - 1.0 : gam;
  return scale * rmean(s, r_floor(x_min, gam, kinetic));
}

[[nodiscard]]
constexpr double
  xmean(const double chi_e, const double gam, const bool kinetic, const double x_min)
{ return xmean(shape(chi_e), gam, kinetic, x_min); }

struct Products {
  /// Emitted photon.
  Vec3 k;
  /// Recoiled lepton.
  Vec3 u1;
};

/// Returns the photon above x_min and the recoiled lepton of lepton u at chi_e with
/// shape s.
///
/// Collinear: the photon leaves along u, the lepton keeps its direction and loses x on
/// the mass shell, the inverse of the self-absorption product.
[[nodiscard]]
constexpr Products
  sample(
    const Shape& s,
    const double chi_e,
    const Vec3& u,
    runko::Rng& rng,
    const bool kinetic,
    const double x_min)
{
  const auto gam = sstd::sqrt(1.0 + toolbox::dot(u, u));
  const auto us  = sstd::sqrt(sstd::max(gam * gam - 1.0, tiny));
  const auto om  = u / us;

  // the deviate covers only the cumulative above the floor: no rejection
  const auto P0 = P(s, r_floor(x_min, gam, kinetic));
  const auto r  = r_of(P0 + (1.0 - P0) * rng.uniform(), s, chi_e);

  const auto x    = (kinetic ? gam - 1.0 : gam) * r;
  const auto gam1 = sstd::max(gam - x, 1.0);
  return Products { x * om, om * sstd::sqrt(sstd::max(gam1 * gam1 - 1.0, 0.0)) };
}

/// Returns the same for lepton u in the field b_perp.
[[nodiscard]]
constexpr Products
  sample(
    const Vec3& u,
    const double b_perp,
    runko::Rng& rng,
    const bool kinetic,
    const double x_min)
{
  const auto chi_e = sstd::sqrt(1.0 + toolbox::dot(u, u)) * b_perp;
  return sample(shape(chi_e), chi_e, u, rng, kinetic, x_min);
}

}  // namespace qed::syn
