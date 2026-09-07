// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <functional>

namespace runko {

/// Constructs YeeLattice to each local tile if it is not already constructed.
void ensure_constructed_yee_lattices(std::reference_wrapper<simulation_context>);

tyvi::actions::sexpr_sender set_EBJ(const tyvi::actions::sexpr&);

}  // namespace runko
