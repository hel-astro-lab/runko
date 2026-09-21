// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "mpi.h"
#include "runko/communication_common.h"
#include "runko/emf/common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/hollow_grid.h"
#include "tyvi/actions.h"

#include <format>


namespace emf {

struct comm_buffs {
  static constexpr auto D = 3uz;
  using value_type        = emf::YeeLattice::value_type;

  using hollow_grid_EB = toolbox::hollow_grid<value_type, D, emf::halo_size>;
  using hollow_grid_J  = toolbox::hollow_grid<value_type, D, 2 * emf::halo_size>;

  hollow_grid_EB E, B;
  hollow_grid_J J;
};

template<runko::comm_mode mode>
runko::exclusive_communicator<mode> make_communicator(runko::simulation_context&);

/// Constructs emf::comm_buffs for virtual and boundary tiles.
void ensure_constructed_comm_buffs(runko::simulation_context&);
/// Updates send buffer in boundary tiles.
void update_send_buff(runko::simulation_context&, runko::comm_mode);

tyvi::actions::sexpr_sender comm_external(runko::simulation_context&, runko::comm_mode);
tyvi::actions::sexpr_sender comm_local(runko::simulation_context&, runko::comm_mode);


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

}  // namespace emf
