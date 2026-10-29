// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "mpi.h"
#include "runko/communication_common.h"
#include "runko/simulation_context.h"
#include "runko/tools/magic_enum/magic_enum_utility.hpp"
#include "tyvi/actions.h"

namespace runko {

/// Creates or gets exclusive_communicator<mode>.
///
/// Due to limitation of MPI (*) all of the communicators are eagerly
/// created at the first invoke. Otherwise parallel make_communicator calls
/// could create invalid communicators as MPI collectives (e.g. MPI_Comm_dub)
/// have to be invoked in the same order on all ranks.
///
/// (*) Lack of tagged MPI_Comm_dub similar to MPI_I{send,recv}.
template<runko::comm_mode mode>
runko::exclusive_communicator<mode> make_communicator(runko::simulation_context&);

/// Constructs communication buffers for virtual and boundary tiles.
void ensure_constructed_comm_buffs(runko::simulation_context&, runko::comm_mode);
/// Updates send buffer in boundary tiles.
void update_send_buff(runko::simulation_context&, runko::comm_mode);

tyvi::actions::sexpr_sender comm_external(runko::simulation_context&, runko::comm_mode);


// Implementation:

template<runko::comm_mode mode>
runko::exclusive_communicator<mode>
  make_communicator(runko::simulation_context& sim)
{
  magic_enum::enum_for_each<runko::comm_mode>([&sim](const auto m_constant) {
    static constexpr runko::comm_mode m = m_constant;
    sim.set_config<runko::exclusive_communicator<m>>([](auto&&) {
      MPI_Comm comm;
      const auto ret = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      if(ret != MPI_SUCCESS) {
        throw std::runtime_error {
          std::format("make_communicator<{}>: MPI_Comm_dub failed", m)
        };
      }
      return runko::exclusive_communicator<m> { comm };
    });
  });

  return sim.get_config<runko::exclusive_communicator<mode>>();
}

}  // namespace runko
