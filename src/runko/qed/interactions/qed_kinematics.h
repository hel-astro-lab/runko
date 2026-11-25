// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "runko/tools/vector.h"

namespace qed {

using Vec3 = toolbox::Vec3<double>;

/// Largest number of proposals of a rejection loop; beyond it the last proposal is
/// taken.
///
/// At acceptance > 0.5 per trip the cap is hit with probability below 1e-19.
inline constexpr int trip_cap = 64;

/// Right-handed triad (i, j, k) with k = kvec and ref in the (i, k) plane.
struct Triad {
  Vec3 i, j, k;
};

[[nodiscard]]
constexpr Triad
  perp_basis(const Vec3& kvec, const Vec3& ref)
{
  auto j  = toolbox::cross(kvec, ref);
  auto nj = toolbox::norm(j);

  if(nj < 1e-10) {
    // ref || kvec: use the Cartesian axis with the smallest component of kvec
    auto alt      = Vec3(0.0, 0.0, 0.0);
    const auto ax = sstd::abs(kvec[0]);
    const auto ay = sstd::abs(kvec[1]);
    const auto az = sstd::abs(kvec[2]);
    if(ax <= ay and ax <= az) {
      alt[0] = 1.0;
    } else if(ay <= az) {
      alt[1] = 1.0;
    } else {
      alt[2] = 1.0;
    }
    j  = toolbox::cross(kvec, alt);
    nj = toolbox::norm(j);
  }

  j = j / nj;
  return Triad { toolbox::cross(j, kvec), j, kvec };
}

/// A four-vector as energy and momentum.
struct FourVec {
  double e;
  Vec3 p;
};

/// Returns Lambda(beta) (p0, p), the components in the frame moving with velocity beta;
/// boost(-beta, .) is the inverse.
[[nodiscard]]
constexpr FourVec
  boost(const Vec3& beta, const double p0, const Vec3& p)
{
  const auto b2  = toolbox::dot(beta, beta);
  const auto g   = 1.0 / sstd::sqrt(1.0 - b2);
  const auto u   = g * beta;
  const auto bp  = toolbox::dot(beta, p);
  const auto fac = b2 > 1e-24 ? (g - 1.0) / b2 : 0.5;  // (g-1)/b^2 -> 1/2 as b -> 0

  const auto p0n = g * p0 - toolbox::dot(u, p);
  const auto pn  = p - p0 * u + (fac * bp) * beta;
  return FourVec { p0n, pn };
}

/// Returns the CM direction of the drawn particle of a pair channel in the triad t.
///
/// The polar axis is t.k, the azimuth origin the plane of t.k and the fixed direction
/// at cosine cfix to it. ntry receives the number of proposals.
[[nodiscard]]
constexpr Vec3
  cm_direction(
    const double cfix,
    const double gcm,
    const double vcm,
    const Triad& t,
    runko::Rng& rng,
    int& ntry)
{
  // The angular factor of e- e+ <-> gamma gamma is F = (gcm^2 + 1) S - 1 - S^2/2 with
  // S = 2/(gcm^2 (1 - vcm^2 w^2)), w the cosine to the fixed direction. Propose w from
  // p(w) ~ 1/(1 - vcm^2 w^2), w = tanh[(2u-1) artanh vcm]/vcm, and accept with F/(C S),
  // C = max F/S over S in [2/gcm^2, 2]; acceptance 0.83 to 1 at every gcm.
  const auto g2  = gcm * gcm;
  const auto fov = [g2](const double S) { return (g2 + 1.0) - 1.0 / S - S / 2.0; };

  // C from the ends and, if reachable, the peak of F/S at S = sqrt 2
  auto C        = sstd::max(fov(2.0 / g2), fov(2.0));
  const auto s2 = 1.4142135623730951;
  if(2.0 / g2 <= s2) { C = sstd::max(C, fov(s2)); }

  const auto at   = sstd::atanh(sstd::min(vcm, 1.0 - 1e-16));
  const auto slow = vcm < 1e-8;  // pair at rest in the CM: w flat
  auto w          = 0.0;
  ntry            = 0;

  while(ntry < trip_cap) {
    ++ntry;
    const auto u1 = rng.uniform();
    w             = slow ? 2.0 * u1 - 1.0
                         : sstd::tanh((2.0 * u1 - 1.0) * at) / sstd::max(vcm, 1e-300);
    const auto S  = 2.0 / (g2 * (1.0 - vcm * vcm * w * w));
    if(fov(S) / C > rng.uniform()) { break; }
  }

  // F depends on w only: the azimuth psi about the fixed direction omf is uniform
  const auto sf  = sstd::sqrt(sstd::max(1.0 - cfix * cfix, 0.0));
  const auto omf = sf * t.i + cfix * t.k;
  const auto e1  = (-cfix) * t.i + sf * t.k;
  const auto psi = sstd::TWOPI * rng.uniform();
  const auto sw  = sstd::sqrt(sstd::max(1.0 - w * w, 0.0));
  return w * omf + (sw * sstd::cos(psi)) * e1 + (sw * sstd::sin(psi)) * t.j;
}

/// Returns an isotropic unit vector.
[[nodiscard]]
constexpr Vec3
  iso_dir(runko::Rng& rng)
{
  const auto cz = -1.0 + 2.0 * rng.uniform();
  const auto ph = sstd::TWOPI * rng.uniform();
  const auto sz = sstd::sqrt(1.0 - cz * cz);
  return Vec3(sz * sstd::cos(ph), sz * sstd::sin(ph), cz);
}

}  // namespace qed
