// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The field the QED channels see. A particle's quantum parameter is chi = (gamma or x)
// b_eff with b_eff = B_eff/B_Q: for the local field the invariant
// sqrt(|E + beta x B|^2 - (beta.E)^2) of the interpolated E, B at the particle, scaled
// from code units by bscale = b_0/B_0; in the field-free box a single uniform value.
// Self-absorption, whose envelope needs one field per tile, takes the tile's median |B|
// across the lepton's pitch angle to the tile's mean field (within 10-20 % of the
// local-field opacity in turbulence).

#include "runko/emf/yee_lattice.h"
#include "runko/qed/interactions/qed_kinematics.h"
#include "runko/tools/math.h"
#include "runko/tools/vector.h"

namespace qed {

/// Linear interpolator of the lattice's E, B at a global position; trivially copyable.
using Interp =
  decltype(std::declval<const emf::YeeLattice&>().interpolate_EB_linear_1st(
    std::declval<std::array<emf::YeeLattice::value_type, 3>>()));

/// The field of a tile as the channels see it.
struct [[nodiscard]] FieldView {
  bool uniform;   // no lattice: every particle sees b0
  double b0;      // the uniform quantum parameter (the box's b_perp)
  double bscale;  // quantum parameter per unit of code field, b_0/B_0
  double b_tile;  // tile's median |B|/B_Q (b0 when uniform); bounds b_ssa
  Vec3 bhat;      // direction of the tile's mean B (z when uniform)
  Interp at;

  /// Field, in B_Q, felt by a particle of velocity beta (a photon: its direction).
  [[nodiscard]] constexpr double
    b_eff(const float x, const float y, const float z, const Vec3& beta) const
  {
    if(uniform) { return b0; }

    const auto eb = at(toolbox::Vec3<float>(x, y, z));
    const auto E  = Vec3(
      static_cast<double>(eb.E[0]),
      static_cast<double>(eb.E[1]),
      static_cast<double>(eb.E[2]));
    const auto B = Vec3(
      static_cast<double>(eb.B[0]),
      static_cast<double>(eb.B[1]),
      static_cast<double>(eb.B[2]));

    const auto f  = E + toolbox::cross(beta, B);
    const auto be = toolbox::dot(beta, E);
    return bscale * sstd::sqrt(sstd::max(toolbox::dot(f, f) - be * be, 0.0));
  }

  /// Field a lepton moving along the unit vector n absorbs in: the tile field across
  /// its pitch angle. The uniform box keeps every lepton perpendicular.
  [[nodiscard]] constexpr double b_ssa(const Vec3& n) const
  { return uniform ? b0 : b_tile * toolbox::norm(toolbox::cross(n, bhat)); }
};

}  // namespace qed
