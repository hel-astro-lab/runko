// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/communication_common.h"

#include <tuple>

/// Parses the 0th and 1st arguments as a rf<simulation_context> and runko::comm_mode.
std::tuple<std::reference_wrapper<runko::simulation_context>, runko::comm_mode>
  runko::args_to_sim_n_comm_mode(const tyvi::actions::sexpr& args)
{
  namespace ta        = tyvi::actions;
  const auto arg_list = std::get<ta::cons>(args);
  const auto arg0     = std::get<ta::atom>(arg_list.car());
  const auto arg_tail = std::get<ta::cons>(arg_list.cdr());
  const auto arg1     = std::get<ta::atom>(arg_tail.car());

  const auto sim  = ta::atom_cast<std::reference_wrapper<simulation_context>>(arg0);
  const auto mode = ta::atom_cast<runko::comm_mode>(arg1);

  if(not sim) {
    throw std::runtime_error(
      "0th argument is not atom containing std::referece_warpper<simulation_context>.");
  }
  if(not mode) {
    throw std::runtime_error("1th argument is not atom containing runko::comm_mode.");
  }

  return std::tuple { sim.value(), mode.value() };
}
