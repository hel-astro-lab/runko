// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/comm/local.h"

#include "runko/comm/cartesian_grid.h"
#include "runko/comm/emf.h"
#include "runko/comm/pic.h"
#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "tyvi/execution.h"

#include <array>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <utility>


namespace runko {

tyvi::actions::sexpr_sender
  comm_local(runko::simulation_context& sim, const runko::comm_mode mode)
{
  namespace te = tyvi::exec;
  namespace ta = tyvi::actions;

  auto f = [sim = std::ref(sim), mode] {
    switch(mode) {
      case runko::comm_mode::emf_E:
      case runko::comm_mode::emf_B:
      case runko::comm_mode::emf_J:
      case runko::comm_mode::emf_J_exchange: return emf::comm_local(sim, mode);
      case runko::comm_mode::pic_particle: return pic::comm_local_particles(sim);
      default:
        throw std::runtime_error {
          std::format("comm_local({}): unhandled comm_mode", mode)
        };
    }
  };
  return te::just() | te::let_value(f);
}

}  // namespace runko
