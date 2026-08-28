// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/pic/virtual_tile.h"
#include "runko/tools/config_parser.h"

#include <array>
#include <cstddef>

namespace qed {

/*! \brief Virtual tile specialization of the qed tile: the pic virtual tile, nothing more. */
template<std::size_t D>
class VirtualTile : virtual public pic::VirtualTile<D> {
public:
  static_assert(D == 3);

  explicit VirtualTile(
    const std::array<std::size_t, 3> tile_indices,
    const toolbox::ConfigParser& config) :
    corgi::Tile<D>(),
    emf::VirtualTile<D>(tile_indices, config),
    pic::VirtualTile<D>(tile_indices, config)
  {
  }

  ~VirtualTile() = default;
};

}  // namespace qed
