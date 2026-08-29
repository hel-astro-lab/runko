// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The field the QED channels see. A particle's quantum parameter is chi = (gamma or x) b_eff,
// with b_eff = B_eff/B_Q the field it feels in units of the critical field: for the local
// field the Lorentz invariant sqrt(|E + beta x B|^2 - (beta.E)^2) of the interpolated E, B at
// the particle, scaled from code units by bscale = b_0/B_0 (b_0 the quantum parameter of the
// reference field B_0); for the field-free box a single uniform value. Self-absorption, whose
// envelope needs one field per tile, takes the tile's median |B| along the lepton's pitch
// angle to the tile's mean field (within 10-20 % of the local-field opacity in turbulence).

#include "runko/pic/tile.h"
#include "runko/qed/interactions/qed_kinematics.h"
#include "runko/tools/math.h"
#include "runko/tools/vector.h"

namespace qed {

/// the tile's linear interpolator of E and B at a global position (a trivially copyable closure)
using Interp = decltype(std::declval<const pic::Tile<3>&>().EB_linear_1st_interpolator());

struct [[nodiscard]] FieldView {
  bool uniform;   // no lattice: every particle sees b0
  double b0;      // the uniform quantum parameter (the box's b_perp)
  double bscale;  // quantum parameter per unit of code field, b_0/B_0
  double b_tile;  // the tile's median |B| in units of B_Q (b0 when uniform); bounds b_ssa
  Vec3 bhat;      // direction of the tile's mean B (z when uniform)
  Interp at;      // E, B at a global position

  /// the field a particle of velocity beta (a photon: its direction) feels, in units of B_Q
  [[nodiscard]] constexpr double b_eff(const float x, const float y, const float z, const Vec3& beta) const
  {
    if(uniform) { return b0; }
    const auto eb = at(toolbox::Vec3<float>(x, y, z));
    const auto E  = Vec3(static_cast<double>(eb.E[0]), static_cast<double>(eb.E[1]), static_cast<double>(eb.E[2]));
    const auto B  = Vec3(static_cast<double>(eb.B[0]), static_cast<double>(eb.B[1]), static_cast<double>(eb.B[2]));
    const auto f  = E + toolbox::cross(beta, B);          // the force field in the particle's frame
    const auto be = toolbox::dot(beta, E);
    return bscale * sstd::sqrt(sstd::max(toolbox::dot(f, f) - be * be, 0.0));
  }

  /// the field a lepton moving along the unit vector n absorbs in: the tile field across its
  /// pitch angle; the uniform box keeps every lepton perpendicular
  [[nodiscard]] constexpr double b_ssa(const Vec3& n) const
  {
    return uniform ? b0 : b_tile * toolbox::norm(toolbox::cross(n, bhat));
  }
};

}  // namespace qed
