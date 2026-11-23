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

std::size_t
  runko::species_n_payload(const toolbox::ConfigParser& config, const std::size_t i)
{
  const auto n = config.get<std::ptrdiff_t>(std::format("n_payload{}", i)).value_or(0);
  if(n < 0) {
    throw std::runtime_error { std::format("n_payload{} = {} is negative.", i, n) };
  }
  return static_cast<std::size_t>(n);
}

double
  runko::species_weight(const toolbox::ConfigParser& config, const std::size_t i)
{
  return config.get<double>(std::format("w{}", i)).value_or(1.0);
}
