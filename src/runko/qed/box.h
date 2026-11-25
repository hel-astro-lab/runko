// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Physical normalization of the QED Monte Carlo of a tile. Weight is in Thomson depths
// across the tile, one unit being W_0 = H^2/sigma_T particles, so tau_T is the summed
// lepton weight. With dt the step in light crossings of the tile, the clocks are
//     N_2 = 1/dt                  a pair reacts w_i w_j sigma-hat F/N_2 times per step
//     N_1 = (lamC/H)/(alpha dt)   a lepton reacts r(chi)/N_1 times per step

#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_ssa.h"

namespace qed {

inline constexpr double LAMC_CM = 3.8615926796e-11;  // reduced Compton wavelength, cm

struct [[nodiscard]] Box {
  double H_lamC;      // tile edge in Compton wavelengths, H/lamC
  double H_box_lamC;  // box edge, which defines ell_B (the tile edge in a 0-D box)
  double dt;          // step in light crossings of the tile, c Delta t/H
  double w_lep;       // weight of every lepton LP, in Thomson depths
  double ell_B;       // magnetic compactness sigma_T H U_B/(m_e c^2); 0: no field
  double x_min;       // floor of the synchrotron band: emission and absorption above it

  [[nodiscard]] constexpr double N2() const { return 1.0 / dt; }
  [[nodiscard]] constexpr double N1() const { return 1.0 / (ALPHA * dt * H_lamC); }

  /// Field B_perp/B_Q from the magnetic compactness of the box.
  [[nodiscard]] constexpr double b_perp() const
  { return ell_B > 0.0 ? b_from_ellB(ell_B, H_box_lamC) : 0.0; }

  /// Energy per step, in W_0 m_e c^2, that the compactness ell injects or removes.
  [[nodiscard]] constexpr double dE(const double ell) const { return ell * dt; }

  /// Escape probability per step of a photon of energy x through a slab of pair depth
  /// tau_T, the step being dt_esc light crossings of the slab.
  [[nodiscard]] static constexpr double
    p_esc(const double x, const double tau_T, const double dt_esc)
  {
    // t_esc = (H/c)(1 + tau_T f(x) sigma_KN); f accounts for the forward scattering of
    // energetic photons and drops to zero above x = 1
    const auto f = x <= 0.1 ? 1.0 : (x < 1.0 ? (1.0 - x) / 0.9 : 0.0);
    return dt_esc / (1.0 + tau_T * f * comp::sigma_kn(x));
  }
};

}  // namespace qed
