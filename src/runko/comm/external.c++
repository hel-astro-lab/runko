// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/comm/external.h"

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/comm/cartesian_grid.h"
#include "runko/comm/emf.h"
#include "runko/comm/pic.h"
#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/vector.h"
#include "tyvi/execution.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <format>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>


namespace runko {

tyvi::actions::sexpr_sender
  comm_external(runko::simulation_context& sim, const runko::comm_mode mode)
{
  switch(mode) {
    case runko::comm_mode::emf_E:
    case runko::comm_mode::emf_B:
    case runko::comm_mode::emf_J: return emf::comm_external(sim, mode);
    case runko::comm_mode::pic_particle: return pic::comm_external_particles(sim);
    default:
      throw std::runtime_error {
        std::format("comm_external({}): unhandled comm_mode", mode)
      };
  }
}

}  // namespace runko
