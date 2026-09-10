// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <functional>

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

/// E -= J for each local tile.
tyvi::actions::sexpr_sender
  add_current(std::reference_wrapper<runko::simulation_context>);

}  // namespace emf
