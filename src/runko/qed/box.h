// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The physical normalization of the QED Monte Carlo of a tile. Weight is measured in units of
// the Thomson depth an LP contributes across the tile, one unit being W_0 = H^2/sigma_T
// particles, so that tau_T is the summed lepton weight and the clocks read
//     N_2 = 1/dt                two-body: a pair reacts w_i w_j sigma-hat F/N_2 times per step,
//     N_1 = (lamC/H)/(alpha dt) one-body: a lepton reacts r(chi)/N_1 times per step,
// with dt the step in light crossings of the tile, and a compactness ell injects (or carries
// out) the energy ell dt per step in units of W_0 m_e c^2.

#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_ssa.h"

namespace qed {

inline constexpr double LAMC_CM = 3.8615926796e-11;  // reduced Compton wavelength, cm

struct [[nodiscard]] Box {
  double H_lamC;  // tile edge in Compton wavelengths, H/lamC
  double dt;      // time step in light crossings of the tile, c Delta t/H
  double w_lep;   // weight of every lepton LP, in Thomson depths
  double ell_B;   // magnetic compactness sigma_T H U_B/(m_e c^2); > 0 opens the field channels
  double x_min;   // floor of the synchrotron band: emission and absorption act above it only

  [[nodiscard]] constexpr double N2() const { return 1.0 / dt; }
  [[nodiscard]] constexpr double N1() const { return 1.0 / (ALPHA * dt * H_lamC); }
  [[nodiscard]] constexpr double b_perp() const
  {
    return ell_B > 0.0 ? b_from_ellB(ell_B, H_lamC) : 0.0;
  }
  [[nodiscard]] constexpr double H_cm() const { return H_lamC * LAMC_CM; }

  /// energy per step, in units of W_0 m_e c^2, that the compactness ell injects
  [[nodiscard]] constexpr double dE(const double ell) const { return ell * dt; }
  /// the compactness a power of dE per step corresponds to
  [[nodiscard]] constexpr double ell_of(const double dE_) const { return dE_ / dt; }

  /// escape probability per step of a photon of energy x through the slab, at pair depth tau_T
  [[nodiscard]] constexpr double p_esc(const double x, const double tau_T) const
  {
    // t_esc = (H/c)(1 + tau_T f(x) sigma_KN), where f accounts for the forward scattering of
    // energetic photons and drops to zero above x = 1
    const auto f = x <= 0.1 ? 1.0 : (x < 1.0 ? (1.0 - x) / 0.9 : 0.0);
    return dt / (1.0 + tau_T * f * comp::sigma_kn(x));
  }
};

}  // namespace qed
