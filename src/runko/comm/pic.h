// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/particles_common.h"
#include "runko/pic/particle.h"
#include "runko/simulation_context.h"
#include "thrust/device_vector.h"
#include "tyvi/actions.h"

#include <map>
#include <span>
#include <vector>

namespace pic {


struct external_comm_buffs {
  using value_type = pic::ParticleContainer::value_type;

  /// All outgoing particle data.
  thrust::device_vector<runko::ParticleState<value_type>> subregion_particle_buff;
  /// One-past-end indices to subregion_particle_buff.
  std::vector<std::size_t> subregion_particle_ends;
};

/// particle type -> incoming particles
using local_comm_buffs = std::map<
  std::size_t,
  std::vector<
    std::span<const runko::ParticleState<pic::external_comm_buffs::value_type>>>>;

tyvi::actions::sexpr_sender comm_external_particles(runko::simulation_context&);
tyvi::actions::sexpr_sender comm_local_particles(runko::simulation_context&);


}  // namespace pic
