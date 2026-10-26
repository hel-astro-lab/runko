// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <optional>
#include <string>

namespace runko {

/// Append `{lap} {<B^2/2>}` to {outdir}/average_B_energy_density.txt.
///
/// Average is taken over all non-halo cells of all tiles.
/// Output directory is read from config (io_outdir), if it is not given explicitly.
tyvi::actions::sexpr_sender average_B_energy_density(
  runko::simulation_context&,
  long lap,
  std::optional<std::string> outdir);

/// Append `{lap} {<E^2/2>}` to {outdir}/average_E_energy_density.txt.
///
/// Average is taken over all non-halo cells of all tiles.
/// Output directory is read from config (io_outdir), if it is not given explicitly.
tyvi::actions::sexpr_sender average_E_energy_density(
  runko::simulation_context&,
  long lap,
  std::optional<std::string> outdir);

/// Append `{lap} {<E_k>_0} {<E_k>_1} ...` to {outdir}/average_kinetic_energy.txt.
///
/// <E_k>_n is the average kinetic energy (in units of m_n c^2) of particles
/// of species n, or zero if there are no such particles.
/// Number of species is read from config (q0/m0, q1/m1, ...).
/// Output directory is read from config (io_outdir), if it is not given explicitly.
tyvi::actions::sexpr_sender average_kinetic_energy(
  runko::simulation_context&,
  long lap,
  std::optional<std::string> outdir);

}  // namespace runko
