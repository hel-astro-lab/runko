// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pybind11/numpy.h"
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

}  // namespace emf
