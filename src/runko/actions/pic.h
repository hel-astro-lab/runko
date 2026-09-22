// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/pic/particle.h"
#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <map>

namespace pic {

using particle_containers = std::map<std::size_t, pic::ParticleContainer>;

/// Constructs ParticleContainers to each local tile if they are not.
void ensure_constructed_particle_containers(runko::simulation_context&);

}  // namespace pic
