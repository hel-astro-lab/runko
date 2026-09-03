// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "tyvi/actions_ast.h"
#include "runko/simulation_context.h"

namespace runko {

enum class symbol : std::uint8_t {
  print,
  println,
  version,
  mt_showcase,
  comm_local,
  comm_external,
  current_context
};

[[nodiscard]]
tyvi::actions::sexpr build_std_env();

[[nodiscard]]
tyvi::actions::sexpr build_sim_env(simulation_context&);

}  // namespace runko
