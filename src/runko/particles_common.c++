// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/particles_common.h"

#include <format>

/// Number of particle species in the config (q0/m0, q1/m1, ...).
std::size_t
  runko::count_species(const toolbox::ConfigParser& config)
{
  auto nspecies = 0uz;
  while(config.contains(std::format("q{}", nspecies)) and
        config.contains(std::format("m{}", nspecies))) {
    ++nspecies;
  }
  return nspecies;
}
