// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/emf/common.h"
#include "runko/emf/yee_lattice.h"
#include "tyvi/actions.h"

namespace emf {

struct comm_buffs {
  static constexpr auto D = 3uz;
  using value_type        = emf::YeeLattice::value_type;

  using hollow_grid_EB = toolbox::hollow_grid<value_type, D, emf::halo_size>;
  using hollow_grid_J  = toolbox::hollow_grid<value_type, D, 2 * emf::halo_size>;

  hollow_grid_EB E, B;
  hollow_grid_J J;
};

tyvi::actions::sexpr_sender comm_external(runko::simulation_context&, runko::comm_mode);
tyvi::actions::sexpr_sender comm_local(runko::simulation_context&, runko::comm_mode);

}  // namespace emf
