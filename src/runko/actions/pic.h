// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pybind11/numpy.h"
#include "runko/particles_common.h"
#include "runko/pic/particle.h"
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
  std::function<runko::ParticleStateBatch(batch_array, batch_array, batch_array)>;

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


/// Push particles in local tiles updating their velocities and positions.
void push_particles(runko::simulation_context&);
/// Deposit current in local tiles from all particls.
void deposit_current(runko::simulation_context&);
/// Sorts the particles in local tiles in order to reduce cache misses.
void sort_particles(runko::simulation_context&);

struct reflectors {
  std::vector<pic::reflector_wall> walls;
};

struct correction_J {
  bool pending;
  using type = runko::VecGrid<emf::YeeLattice::value_type>;
  type J;
};

/// Registers the given reflector to the context variables of type pic::reflectors.
///
/// see: https://github.com/skypjack/entt/wiki/Entity-Component-System#context-variables
void register_reflector_wall(runko::simulation_context&, const pic::reflector_wall&);

/// Reflect particles that crossed any registered reflectors in local tiles.
///
/// Must be called after push_particles and before deposit_current.
/// Modifies particle positions/velocities for reflected particles
/// and stores correction currents that deposit_current will add.
void reflect_particles(runko::simulation_context&);

/// Update reflector wall locations by their velocity * cfl.
void advance_reflector_walls(runko::simulation_context&);

}  // namespace pic
