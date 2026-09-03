// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "entt/entt.hpp"
#include "mpi.h"
#include "runko/communication_common.h"
#include "runko/emf/common.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/vector.h"
#include "tyvi/actions_ast.h"
#include "tyvi/mdspan.h"
#include "tyvi/sstd.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <format>
#include <map>
#include <ranges>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace runko {

/// Marks local tiles.
struct local_tile_tag {};

struct comm_info_type {
  int rank { -1 };
  int tag { -1 };
};

/// Marks boundary tiles.
///
/// Boundary tiles have a corresponding virtual tiles (dest) somewhere.
struct boundary_tile_tag {
  std::vector<comm_info_type> dest_infos;
};

/// Marks virtual tile.
///
/// Virtual tiles have a corresponding boundary tile (source) somewhere.

struct virtual_tile_tag {
  comm_info_type source_info;
};

template<std::size_t rank>
using cartesian_index = toolbox::VecD<std::ptrdiff_t, rank>;

/// View over indices in Moore neighborhood around given index.
template<std::size_t rank>
std::ranges::view auto moore_neighs(const cartesian_index<rank>&);

/// Similar to moore_neighs but adds direction to each index.
///
/// value_type is std::tuple<cartesian_index<rank>, grid_neighbor<rank>>.
template<std::size_t rank>
std::ranges::view auto moore_neighs_n_dirs(const cartesian_index<rank>&);

template<std::size_t rank>
struct cartesian_neighbors {
  void set(const runko::grid_neighbor<rank>&, entt::registry::entity_type);
  auto get(const runko::grid_neighbor<rank>&) const
    -> std::optional<entt::registry::entity_type>;

private:
  std::map<runko::grid_neighbor<rank>, entt::registry::entity_type> neighbors_;
};


struct simulation_context {
  /// FIXME: make all view_tiles usage const correct
  std::ranges::view auto view_tiles();
  template<typename... Ts>
  std::ranges::view auto view_tiles();
  template<typename... Ts>
  std::ranges::view auto view_tiles() const;

  const toolbox::ConfigParser config;
  entt::registry tiles;
  using tile_id_type = entt::registry::entity_type;
};

template<std::size_t rank>
auto wrap_cartesian_index(const cartesian_index<rank>&);

/// Computes and sets cartesian neighbors.
///
/// For each local tile finds and sets its cartesian_neighbors.
/// If a neighboring tile is missing, it will be constructed as a virtual tile.
/// Adds boundary_tile_tags for the boundary tiles.
/// However, boundary_tile_tag::dest_infos and virtual_tile_tag::source_rank
/// are left uninitialized.
template<std::size_t rank>
void set_cartesian_neighbors(simulation_context&);


/// Computes cartesian tile rank info.
///
/// Communicates all virtual tiles with cartesian_index<rank> in MPI_COMM_WORLD
/// and sets communication infos for boundary and virtual tiles.
/// Assumes that all cartesian tiles are located at some rank.
template<std::size_t rank>
void set_cartesian_comm_infos(simulation_context&);

/// Maps cartesian index (potentially from virtual tile) to a MPI tag in MPI_COMM_WORLD.
template<std::size_t rank>
int cartesian_index_to_mpi_tag(
  const std::vector<std::ptrdiff_t>& n_tiles,
  const cartesian_index<rank>&);

auto sim_env(simulation_context&) -> tyvi::actions::sexpr;


// Implementation:

template<std::size_t rank>
void
  cartesian_neighbors<rank>::set(
    const runko::grid_neighbor<rank>& dir,
    const entt::registry::entity_type id)
{ this->neighbors_[dir] = id; }

template<std::size_t rank>
auto
  cartesian_neighbors<rank>::get(const runko::grid_neighbor<rank>& dir) const
  -> std::optional<entt::registry::entity_type>
{
  if(this->neighbors_.contains(dir)) { return this->neighbors_.at(dir); }
  return {};
}

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles()
{ return this->tiles.view<Ts...>().each(); }

template<typename... Ts>
std::ranges::view auto
  simulation_context::view_tiles() const
{ return this->tiles.view<const Ts...>().each(); }

template<std::size_t rank>
void
  set_cartesian_neighbors(simulation_context& sim)
{

  using index_type = runko::cartesian_index<rank>;
  using id_type    = decltype(sim.tiles)::entity_type;

  auto index_to_id = std::map<index_type, id_type> {};
  // Virtual tiles to be constructed.
  auto required = std::set<index_type> {};

  namespace rn = std::ranges;
  namespace rv = std::views;

  const auto n_tiles = toolbox::get_extent_list(sim.config, "n_tiles", rank);

  for(const auto& [id, index]: sim.view_tiles<index_type>()) {
    index_to_id[index] = id;
  }

  for(const auto& [_, index]: sim.view_tiles<index_type, local_tile_tag>()) {
    for(const auto& x: moore_neighs(index)) { std::ignore = required.insert(x); }
  }

  auto to_be_constructed =
    required | rv::filter([&](const auto& x) { return not index_to_id.contains(x); });
  for(const auto& x: to_be_constructed) {
    const auto id  = sim.tiles.create();
    index_to_id[x] = id;
    sim.tiles.emplace<virtual_tile_tag>(id);
    sim.tiles.emplace<index_type>(id, x);
  }

  for(const auto& [id, index]: sim.view_tiles<index_type, local_tile_tag>()) {
    auto& neighbors =
      sim.tiles.emplace_or_replace<runko::cartesian_neighbors<rank>>(id);
    for(const auto& [x, dir]: moore_neighs_n_dirs(index)) {
      const auto neigh_id = index_to_id.at(x);
      neighbors.set(dir, neigh_id);

      const auto neigh_is_virtual =
        nullptr != sim.tiles.try_get<virtual_tile_tag>(neigh_id);
      if(neigh_is_virtual) {
        std::ignore = sim.tiles.emplace_or_replace<boundary_tile_tag>(id);
      }
    }
  }
}
template<std::size_t rank>
auto
  wrap_cartesian_index(
    const cartesian_index<rank>& idx,
    const std::vector<std::ptrdiff_t>& n_tiles)
{

  return [&]<std::size_t... I>(std::index_sequence<I...>) {
    // Index can be negative, so we can not just use remainder operator (%),
    // as we want to calculate proper modulo.
    return cartesian_index<rank>(
      (((idx[I] % n_tiles[I]) + n_tiles[I]) % n_tiles[I])...);
  }(std::make_index_sequence<rank>());
}

template<std::size_t rank>
void
  set_cartesian_comm_infos(simulation_context& sim)
{

  using index_type = runko::cartesian_index<rank>;

  namespace rn = std::ranges;
  namespace rv = std::views;


  const auto n_tiles = toolbox::get_extent_list(sim.config, "n_tiles", rank);

  auto my_local_indices = std::vector<index_type> {};
  auto my_virt_indices  = std::vector<index_type> {};

  for(const auto& [_, index]: sim.view_tiles<index_type, local_tile_tag>()) {
    my_local_indices.push_back(index);
  }

  for(const auto& [_, index, __]: sim.view_tiles<index_type, virtual_tile_tag>()) {
    my_virt_indices.push_back(index);
  }

  // FIXME: senderify this
  /* We want to gather local_offsets from every rank,
     but first we have to figure out how many tiles each rank has. */

  int comm_rank, comm_size;

  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank)) {
    throw std::runtime_error { "MPI_Comm_rank(MPI_COMM_WORLD, ...) failed." };
  }
  if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &comm_size)) {
    throw std::runtime_error { "MPI_Comm_size(MPI_COMM_WORLD, ...) failed." };
  }

  // This is a bit pendantic :)
  static_assert(sizeof(std::size_t) == 8uz);
  static_assert(sizeof(int) == 4uz);
  static_assert(CHAR_BIT == 8uz);

  auto local_tile_counts = std::vector<int>(static_cast<std::size_t>(comm_size));
  const auto my_local_tile_count = checked_cast<int>(my_local_indices.size());

  if(
    MPI_SUCCESS != MPI_Allgather(
                     reinterpret_cast<const void*>(&my_local_tile_count),
                     1,
                     MPI_INT32_T,
                     reinterpret_cast<void*>(local_tile_counts.data()),
                     1,
                     MPI_INT32_T,
                     MPI_COMM_WORLD)) {
    throw std::runtime_error { "MPI_Allgather failed (local indices)." };
  }

  auto virt_tile_counts         = std::vector<int>(static_cast<std::size_t>(comm_size));
  const auto my_virt_tile_count = checked_cast<int>(my_virt_indices.size());

  if(
    MPI_SUCCESS != MPI_Allgather(
                     reinterpret_cast<const void*>(&my_virt_tile_count),
                     1,
                     MPI_INT32_T,
                     reinterpret_cast<void*>(virt_tile_counts.data()),
                     1,
                     MPI_INT32_T,
                     MPI_COMM_WORLD)) {
    throw std::runtime_error { "MPI_Allgather failed (virtual indices)." };
  }

  auto local_displs = std::vector<int>(static_cast<std::size_t>(comm_size));
  auto virt_displs  = std::vector<int>(static_cast<std::size_t>(comm_size));

  local_displs.front() = 0;
  for(auto i = 1uz; i < local_displs.size(); ++i) {
    local_displs[i] = local_displs[i - 1uz] + local_displs[i];
  }
  virt_displs.front() = 0;
  for(auto i = 1uz; i < virt_displs.size(); ++i) {
    virt_displs[i] = virt_displs[i - 1uz] + virt_displs[i];
  }

  auto local_indices = std::vector<index_type>(
    static_cast<std::size_t>(local_displs.back() + local_tile_counts.back()));
  auto virt_indices = std::vector<index_type>(
    static_cast<std::size_t>(virt_displs.back() + virt_tile_counts.back()));

  /* We don't create a custom MPI datatype for index_type.
     We'll just handle bytes, so this has to be taken into account,
     in recvcount and displs. */

  for(auto& x: local_tile_counts) { x *= static_cast<int>(sizeof(index_type)); };
  for(auto& x: virt_tile_counts) { x *= static_cast<int>(sizeof(index_type)); };


  if(
    const auto err = MPI_Allgatherv(
      reinterpret_cast<const void*>(my_local_indices.data()),
      checked_cast<int>(sizeof(index_type) * my_local_indices.size()),
      MPI_BYTE,
      reinterpret_cast<void*>(local_indices.data()),
      local_tile_counts.data(),
      local_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
    err != MPI_SUCCESS) {
    throw std::runtime_error {
      std::format("MPI_Allgatherv failed (local indices): {}", err)
    };
  }

  if(
    const auto err = MPI_Allgatherv(
      reinterpret_cast<const void*>(my_virt_indices.data()),
      checked_cast<int>(sizeof(index_type) * my_virt_indices.size()),
      MPI_BYTE,
      reinterpret_cast<void*>(virt_indices.data()),
      virt_tile_counts.data(),
      virt_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
    err != MPI_SUCCESS) {
    throw std::runtime_error {
      std::format("MPI_Allgatherv failed (virt indices): {}", err)
    };
  }

  auto local_index_to_rank = std::map<index_type, int> {};
  {
    auto current_rank = 0;
    for(const auto& [n, idx]: std::views::enumerate(local_indices)) {
      if(
        current_rank + 1 < std::ranges::ssize(local_displs) and
        n >= local_displs.at(static_cast<std::size_t>(current_rank + 1))) {
        current_rank += 1;
      }
      local_index_to_rank[idx] = current_rank;
    }
  }
  auto virt_index_to_rank = std::map<index_type, int> {};
  {
    auto current_rank = 0;
    for(const auto& [n, idx]: std::views::enumerate(virt_indices)) {
      if(
        current_rank + 1 < std::ranges::ssize(virt_displs) and
        n >= local_displs.at(static_cast<std::size_t>(current_rank + 1))) {
        current_rank += 1;
      }
      virt_index_to_rank[idx] = current_rank;
    }
  }

  for(auto&& [_, idx, virt]: sim.view_tiles<index_type, virtual_tile_tag>()) {
    virt.source_info = comm_info_type {
      .rank = local_index_to_rank.at(wrap_cartesian_index<rank>(idx, n_tiles)),
      .tag  = cartesian_index_to_mpi_tag<rank>(n_tiles, idx)
    };
  }

  for(auto&& [_, idx, boundary]: sim.view_tiles<index_type, boundary_tile_tag>()) {
    boundary.dest_infos.clear();

    for(const auto& [virt_idx, virt_rank]: virt_index_to_rank) {
      if(idx == wrap_cartesian_index<rank>(virt_idx, n_tiles)) {
        boundary.dest_infos.push_back(
          comm_info_type { .rank = virt_rank,
                           .tag =
                             cartesian_index_to_mpi_tag<rank>(n_tiles, virt_idx) });
      }
    }
  }
}

template<std::size_t rank>
int
  cartesian_index_to_mpi_tag(
    const std::vector<std::ptrdiff_t>& n_tiles,
    const cartesian_index<rank>& idx)
{
  if(n_tiles.size() != rank) {
    throw std::runtime_error { "cartesian_index_to_mpi_tag: rank != n_tiles.size()" };
  }

  using M = std::layout_left::mapping<std::dextents<std::size_t, rank>>;

  const auto tag = runko::checked_cast<int>([&]<std::size_t... I>(
                                              std::index_sequence<I...>) {
    /// cartesian indices of a virtual tiles are in range [-1, N],
    /// while indices of a local tiles are in range [0, N).
    /// In order to support virtual tiles, we have to shift the indices by one
    /// and widen the extents by two.
    return M(std::dextents<std::ptrdiff_t, rank>((n_tiles[I] + 2)...))((idx[I] + 1)...);
  }(std::make_index_sequence<rank>()));

  static int mpi_tag_ub = -1;
  if(mpi_tag_ub == -1) {
    void* loc = nullptr;
    int flag  = false;
    if(MPI_SUCCESS != MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &loc, &flag)) {
      throw std::runtime_error {
        "MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, ...) failed."
      };
    }

    if(not flag) {
      throw std::logic_error {
        "MPI_TAG_UB is defined for MPI_COMM_WORLD but did not have value."
      };
    }

    mpi_tag_ub = *reinterpret_cast<int*>(loc);

    if(mpi_tag_ub == -1) {
      throw std::logic_error {
        "MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, ...) sanity check failed."
      };
    }
  }

  if(tag < 0 or tag >= mpi_tag_ub) {
    throw std::runtime_error(
      std::format(
        "Calculated tag outside of [0, MPI_TAG_UB = {}): {}",
        mpi_tag_ub,
        tag));
  }

  return tag;
}

template<std::size_t rank>
std::ranges::view auto
  moore_neighs(const cartesian_index<rank>& idx)
{
  namespace rv     = std::views;
  using neigh_type = runko::grid_neighbor<rank>;

  return rv::iota(0uz, std::pow(3uz, rank)) | rv::transform([idx](const auto n) {
           return idx + neigh_type::from_index(n).template to_vec<std::ptrdiff_t>();
         }) |
         rv::filter([idx](const auto& x) { return x != idx; });
}

template<std::size_t rank>
std::ranges::view auto
  moore_neighs_n_dirs(const cartesian_index<rank>& idx)
{
  namespace rv     = std::views;
  using neigh_type = runko::grid_neighbor<rank>;

  return rv::iota(0uz, std::pow(3uz, rank)) | rv::transform([idx](const auto n) {
           const auto dir  = neigh_type::from_index(n);
           const auto vdir = dir.template to_vec<std::ptrdiff_t>();
           return std::tuple { idx + vdir, dir };
         }) |
         rv::filter(
           [](const auto& x) { return std::get<1>(x) != grid_neighbor_origo<rank>; });
}

}  // namespace runko
