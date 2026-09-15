// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pybind11/numpy.h"
#include "runko/emf/antenna.h"
#include "runko/emf/edge_bc.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <functional>
#include <vector>

namespace emf {

/// Constructs YeeLattice to each local tile if it is not already constructed.
void ensure_constructed_yee_lattices(std::reference_wrapper<runko::simulation_context>);

using vector_field_function =
  std::function<std::tuple<double, double, double>(double, double, double)>;

/// Sets E, B and J for each local tile.
tyvi::actions::sexpr_sender set_EBJ(
  std::reference_wrapper<runko::simulation_context>,
  vector_field_function E,
  vector_field_function B,
  vector_field_function J);

using batch_array = pybind11::array_t<double>;
using batch_vector_field_function =
  std::function<batch_array(batch_array, batch_array, batch_array)>;

/// Batch sets E, B and J for each local tile.
tyvi::actions::sexpr_sender batch_set_EBJ(
  std::reference_wrapper<runko::simulation_context>,
  batch_vector_field_function Ex,
  batch_vector_field_function Ey,
  batch_vector_field_function Ez,
  batch_vector_field_function Bx,
  batch_vector_field_function By,
  batch_vector_field_function Bz,
  batch_vector_field_function Jx,
  batch_vector_field_function Jy,
  batch_vector_field_function Jz);

/// E -= J for each local tile.
tyvi::actions::sexpr_sender
  add_current(std::reference_wrapper<runko::simulation_context>);

struct antennas {
  std::vector<emf::antenna_mode> modes;
};

/// Registers a given antenna to as a entt::registry context variable.
///
/// see: https://github.com/skypjack/entt/wiki/Entity-Component-System#context-variables
void register_antenna(runko::simulation_context&, emf::antenna_mode);

struct antenna_buffers {
  runko::VecGrid<emf::YeeLattice::value_type> vec_pot;
  runko::VecGrid<emf::YeeLattice::value_type> generated_B;
};

/// Deposits current from registered antennas to local tiles.
///
/// If emf::antenna_mode contains modes with lap_coeffs,
/// uses the latest one and remove it.
tyvi::actions::sexpr_sender deposit_antenna_current(runko::simulation_context&);

tyvi::actions::sexpr_sender push_e(runko::simulation_context&);
tyvi::actions::sexpr_sender push_half_b(runko::simulation_context&);
tyvi::actions::sexpr_sender filter_current(runko::simulation_context&);

struct boundary_conditions {
  std::vector<emf::edge_bc> edges;
};

/// Registers a given edge_bc to as a entt::registry context variable.
///
/// Type of the context variable is emf::boundary_conditions.
/// see: https://github.com/skypjack/entt/wiki/Entity-Component-System#context-variables
void register_edge_bc(runko::simulation_context&, const emf::edge_bc&);

/// Applies given edge boundary condition to local tiles.
tyvi::actions::sexpr_sender
  apply_edge_bc(runko::simulation_context&, const emf::edge_bc&, runko::comm_mode);

/// Applies registered edge boundary conditions to local tiles.
tyvi::actions::sexpr_sender
  apply_edge_bcs(runko::simulation_context&, runko::comm_mode);


}  // namespace emf
