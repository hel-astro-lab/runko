// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Frame machinery shared by the two-body QED samplers: the Lorentz boost, the triad the
// differential cross sections are written in, and the CM angular draw of the pair channels.
// None of it knows anything about the Monte Carlo scheme. Units: m_e c^2, m_e c.

#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "runko/tools/vector.h"

namespace qed {

using Vec3 = toolbox::Vec3<double>;

/// Largest number of proposals of a rejection loop; beyond it the last proposal is taken.
/// At the measured acceptance rates (> 0.5 per trip) the cap is reached with probability
/// below 1e-19, so it only bounds the work of a lane.
inline constexpr int TRIP_CAP = 64;

/// Right-handed triad (i, j, k) with k = kvec and ref in the (i, k) plane; handles ref || kvec.
struct [[nodiscard]] Triad {
  Vec3 i, j, k;
};

[[nodiscard]]
constexpr Triad
  perp_basis(const Vec3& kvec, const Vec3& ref)
{
  // kvec is the polar axis of the angular sampler and ref fixes the azimuth origin, so
  // that phi is measured from the plane the two of them span.
  auto j  = toolbox::cross(kvec, ref);  // normal of the (kvec, ref) plane
  auto nj = toolbox::norm(j);
  if(nj < 1e-10) {
    // ref parallel to kvec leaves the plane undefined; any perpendicular direction serves,
    // and the Cartesian axis kvec leans on least is guaranteed not to be parallel to it.
    // pick the Cartesian axis with the smallest component of kvec
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
  j = j / nj;                            // unit normal of the plane
  return Triad { toolbox::cross(j, kvec), j, kvec };  // i completes the right-handed triad
}

/// A four-vector as energy and momentum.
struct [[nodiscard]] FourVec {
  double e;
  Vec3 p;
};

/// Lambda(beta) applied to (p0, p): the components in the frame moving with velocity beta.
/// boost(-beta, .) is the inverse transformation.
[[nodiscard]]
constexpr FourVec
  boost(const Vec3& beta, const double p0, const Vec3& p)
{
  // ---- boost quantities
  const auto b2  = toolbox::dot(beta, beta);           // beta^2
  const auto g   = 1.0 / sstd::sqrt(1.0 - b2);         // Lorentz factor of the boost
  const auto u   = g * beta;                           // four-velocity gamma beta
  const auto bp  = toolbox::dot(beta, p);              // projection of p on the boost axis
  const auto fac = b2 > 1e-24 ? (g - 1.0) / b2 : 0.5;  // (g-1)/b^2 -> 1/2 as beta -> 0
  const auto p0n = g * p0 - toolbox::dot(u, p);        // transformed energy
  const auto pn  = p - p0 * u + (fac * bp) * beta;     // only the longitudinal part moves
  return FourVec { p0n, pn };
}

/// CM direction of the drawn particle of a conversion channel, in the triad whose polar axis
/// is t.k and whose azimuth origin is the plane of t.k and the fixed direction, the latter at
/// cosine cfix to the axis. ntry receives the number of proposals.
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
  // The angular factor of e- e+ <-> gamma gamma depends on the CM angle only through
  //     S = 2/(gcm^2 (1 - vcm^2 w^2)),  w = cos(angle to the fixed direction),
  // as F = (gcm^2 + 1) S - 1 - S^2/2. Proposing w from S itself -- p(w) ~ 1/(1 - vcm^2 w^2),
  // inverse cumulative w = tanh[(2u-1) artanh vcm]/vcm -- and thinning with F/(C S), C the
  // maximum of F/S over the reachable S in [2/gcm^2, 2], accepts 0.83 to 1 at every gcm.
  // ---- the bound C of F/S over the reachable S: its ends and, if inside, its peak
  const auto g2  = gcm * gcm;                         // gamma_cm^2
  const auto fov = [g2](const double S) { return (g2 + 1.0) - 1.0 / S - S / 2.0; };  // F/S
  auto C         = sstd::max(fov(2.0 / g2), fov(2.0));
  const auto s2  = 1.4142135623730951;            // F/S peaks at S = sqrt 2
  if(2.0 / g2 <= s2) { C = sstd::max(C, fov(s2)); }
  const auto at   = sstd::atanh(sstd::min(vcm, 1.0 - 1e-16));  // artanh(beta_cm) of the proposal
  const auto slow = vcm < 1e-8;                   // pair at rest in the CM: flat in w
  auto w          = 0.0;
  ntry            = 0;
  // rejection loop: propose w from 1/(1 - vcm^2 w^2), accept with F/(C S)
  while(ntry < TRIP_CAP) {
    ++ntry;
    const auto u1 = rng.uniform();
    w = slow ? 2.0 * u1 - 1.0 : sstd::tanh((2.0 * u1 - 1.0) * at) / sstd::max(vcm, 1e-300);
    const auto S = 2.0 / (g2 * (1.0 - vcm * vcm * w * w));  // the invariant of the proposal
    if(fov(S) / C > rng.uniform()) { break; }    // = F/(C S), in [0, 1]
  }
  // w is the cosine to the fixed direction, whose azimuth about it is uniform because F
  // depends on nothing else; (e1, t.j, omf) is an orthonormal triad built on that direction.
  const auto sf  = sstd::sqrt(sstd::max(1.0 - cfix * cfix, 0.0));  // sine to the fixed direction
  const auto omf = sf * t.i + cfix * t.k;         // the fixed direction in the triad
  const auto e1  = (-cfix) * t.i + sf * t.k;      // perpendicular to it, in the (i, k) plane
  const auto psi = 6.283185307179586 * rng.uniform();  // azimuth about the fixed direction
  const auto sw  = sstd::sqrt(sstd::max(1.0 - w * w, 0.0));  // sine of the drawn angle
  return w * omf + (sw * sstd::cos(psi)) * e1 + (sw * sstd::sin(psi)) * t.j;
}

/// Isotropic unit vector.
[[nodiscard]]
constexpr Vec3
  iso_dir(runko::Rng& rng)
{
  const auto cz = -1.0 + 2.0 * rng.uniform();     // cos(theta), uniform on the sphere
  const auto ph = 6.283185307179586 * rng.uniform();  // azimuth
  const auto sz = sstd::sqrt(1.0 - cz * cz);      // sin(theta)
  return Vec3(sz * sstd::cos(ph), sz * sstd::sin(ph), cz);
}

}  // namespace qed
