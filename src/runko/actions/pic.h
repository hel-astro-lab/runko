// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pybind11/numpy.h"
#include "runko/particles_common.h"
#include "runko/pic/particle.h"
#include "runko/pic/tile.h"
#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <functional>
#include <map>

namespace pic {

using particle_containers = std::map<std::size_t, pic::ParticleContainer>;

/// Constructs ParticleContainers to each local tile if they are not.
void ensure_constructed_particle_containers(runko::simulation_context&);

using particle_generator =
  std::function<std::vector<runko::ParticleState<double>>(double, double, double)>;

/// Inject particles to local tiles based on given generator.
///
/// Generator is called for each cell coordinate.
///
/// Particle type is assumed to be configured.
void inject_to_each_cell(
  runko::simulation_context&,
  std::size_t particle_type,
  particle_generator);

/// Inject copy of given particles all local tiles.
///
/// Should be used with caution as this does not really make sense
/// if there are more than one local tile.
///
/// Particle type is assumed to be configured.
void inject(
  runko::simulation_context&,
  std::size_t particle_type,
  std::vector<runko::ParticleState<double>>);

using batch_array = pybind11::array_t<double>;
using batch_particle_generator =
  std::function<pic::ParticleStateBatch(batch_array, batch_array, batch_array)>;

/// Inject particles to local tiles based on given generator.
///
/// Generator is called once with all cell coordinates.
///
/// Particle type is assumed to be configured.
void batch_inject_to_cells(
  runko::simulation_context&,
  std::size_t particle_type,
  batch_particle_generator);

/// Inject particles in a stripe between x_left and x_right to local tiles.
///
/// Cells overlapping [x_left, x_right) are passed to the generator; generated
/// particles outside [x_left, x_right) are dropped, so partial edge cells get
/// the matching fraction. Full y and z extent of the tile is used.
void batch_inject_in_x_stripe(
  runko::simulation_context&,
  std::size_t particle_type,
  batch_particle_generator pgen,
  double x_left,
  double x_right);

}  // namespace pic
