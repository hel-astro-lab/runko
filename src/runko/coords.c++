// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/coords.h"

#include "runko/tools/config_parser.h"

namespace runko {

global_coordinates_closure<3>
  global_coordinates(const simulation_context& sim, const std::array<double, 3>& idx)
{
  const auto cells  = toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);
  const auto vcells = toolbox::Vec3<double>(cells[0], cells[1], cells[2]);
  const auto vidx   = toolbox::Vec3<double>(idx);

  return { vcells, vidx };
}

}  // namespace runko
