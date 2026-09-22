// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "mpi.h"
#include "runko/communication_common.h"
#include "runko/simulation_context.h"
#include "tyvi/actions.h"

namespace runko {

template<runko::comm_mode mode>
runko::exclusive_communicator<mode> make_communicator(runko::simulation_context&);

/// Constructs communication buffers for virtual and boundary tiles.
void ensure_constructed_comm_buffs(runko::simulation_context&);
/// Updates send buffer in boundary tiles.
void update_send_buff(runko::simulation_context&, runko::comm_mode);

tyvi::actions::sexpr_sender comm_external(runko::simulation_context&, runko::comm_mode);


// Implementation:

template<runko::comm_mode mode>
runko::exclusive_communicator<mode>
  make_communicator(runko::simulation_context& sim)
{
  return sim.get_n_set_config<runko::exclusive_communicator<mode>>([](auto&&) {
    MPI_Comm comm;
    const auto ret = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    if(ret != MPI_SUCCESS) {
      throw std::runtime_error {
        std::format("make_communicator<{}>: MPI_Comm_dub failed", mode)
      };
    }
    return runko::exclusive_communicator<mode> { comm };
  });
}

}  // namespace runko
