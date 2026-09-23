// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

namespace runko {

enum class symbol : std::uint8_t {
  print,
  println,
  format,
  version,
  mt_showcase,
  sequence,
  comm_local,
  comm_external,
  current_context,
  ensure_constructed_yee_lattices,
  set_EBJ,
  batch_set_EBJ,
  add_current,
  register_antenna,
  deposit_antenna_current,
  push_e,
  push_half_b,
  filter_current,
  register_edge_bc,
  apply_edge_bc,
  apply_edge_bcs,
  clear_bcs,
  ensure_constructed_particle_containers,
  inject_to_each_cell,
  inject,
  batch_inject_to_cells,
  batch_inject_in_x_stripe,
  push_particles,
  set_cartesian_neighbors,
  set_cartesian_comm_infos
};

[[nodiscard]]
tyvi::actions::sexpr build_std_env();

[[nodiscard]]
tyvi::actions::sexpr build_sim_env(simulation_context&);

}  // namespace runko
