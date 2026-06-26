// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/emf/virtual_tile.h"

#include "runko/communication_common.h"
#include "runko/tools/system.h"
#include "thrust/host_vector.h"

namespace emf {

template<std::size_t D>
VirtualTile<D>::VirtualTile(
  const std::array<std::size_t, 3> tile_indices,
  const toolbox::ConfigParser& p) :
  corgi::Tile<D>()
{

  const auto cells = p.get_or_throw<std::vector<std::ptrdiff_t>>("n_cells_per_tile");
  this->extents_wout_halo_ = { static_cast<std::size_t>(cells[0]),
                               static_cast<std::size_t>(cells[1]),
                               static_cast<std::size_t>(cells[2]) };

  E_ = hollow_grid_EB(extents_wout_halo());
  B_ = hollow_grid_EB(extents_wout_halo());
  J_ = hollow_grid_J(extents_with_halo());

  const auto xmin = 0.0;
  const auto ymin = 0.0;
  const auto zmin = 0.0;

  const auto [i, j, k] = tile_indices;

  const auto tiles = p.get_or_throw<std::vector<std::ptrdiff_t>>("n_tiles");
  const auto Nx    = static_cast<std::size_t>(tiles[0]);
  const auto Ny    = static_cast<std::size_t>(tiles[1]);
  const auto Nz    = static_cast<std::size_t>(tiles[2]);

  if(Nx <= i or Ny <= j or Nz <= k) {
    throw std::runtime_error { "Trying to create tile outside of configured grid." };
  }

  this->index = { tile_indices[0], tile_indices[1], tile_indices[2] };

  // Tile size.
  const auto [Lx, Ly, Lz] = extents_wout_halo_;

  this->set_tile_mins(
    { xmin + static_cast<double>(i * Lx),
      ymin + static_cast<double>(j * Ly),
      zmin + static_cast<double>(k * Lz) });
  this->set_tile_maxs(
    { xmin + static_cast<double>((i + 1uz) * Lx),
      ymin + static_cast<double>((j + 1uz) * Ly),
      zmin + static_cast<double>((k + 1uz) * Lz) });
}

template<std::size_t D>
std::array<std::size_t, 3>
  VirtualTile<D>::extents_wout_halo() const
{
  return this->extents_wout_halo_;
}

template<std::size_t D>
std::array<std::size_t, 3>
  VirtualTile<D>::extents_with_halo() const
{

  return { this->extents_wout_halo_[0] + 2 * halo_size,
           this->extents_wout_halo_[1] + 2 * halo_size,
           this->extents_wout_halo_[2] + 2 * halo_size };
}

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  VirtualTile<D>::recv_data(
    mpi4cpp::mpi::communicator& comm,
    const int orig,
    const int mode,
    const int tag)
{

#ifndef TYVI_BACKEND_CPU
  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error { "GPU backend requires GPU-aware MPI." };
  }
#endif

  using runko::comm_mode;

  auto make_irecv = [&](const auto s) {
    return comm.irecv(orig, tag, s.data(), runko::checked_cast<int>(s.size()));
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E: return { make_irecv(E_.span()) };
    case comm_mode::emf_B: return { make_irecv(B_.span()) };
    case comm_mode::emf_J: return { make_irecv(J_.span()) };
    default:
      throw std::logic_error { std::format(
        "emf::Tile::recv_data does not support given communication mode: {}",
        mode) };
  }
}
}  // namespace emf

template class emf::VirtualTile<3>;
