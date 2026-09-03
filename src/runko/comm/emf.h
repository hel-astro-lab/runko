// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pika/execution.hpp"
#include "pika/mpi.hpp"
#include "runko/emf/common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/hollow_grid.h"
#include "runko/tools/vector.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <format>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace runko {

struct emf_comm_buffs {
  static constexpr auto D = 3uz;
  using value_type        = emf::YeeLattice::value_type;

  using hollow_grid_EB = toolbox::hollow_grid<value_type, D, emf::halo_size>;
  using hollow_grid_J  = toolbox::hollow_grid<value_type, D, 2 * emf::halo_size>;

  hollow_grid_EB E, B;
  hollow_grid_J J;
};

/// Constructs emf_comm_buffs for virtual and boundary tiles.
auto ensure_constructed_emf_comm_buffs(std::reference_wrapper<simulation_context>);
/// Updates B send buffer in boundary tiles.
auto update_send_buff_B(std::reference_wrapper<simulation_context>);
auto comm_external_B(std::reference_wrapper<simulation_context>);
auto comm_local_B(std::reference_wrapper<simulation_context>);

// Implementation:

auto
  ensure_constructed_emf_comm_buffs(
    const std::reference_wrapper<simulation_context> sim)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  const auto n_cells =
    toolbox::get_extent_list(sim.get().config, "n_cells_per_tile", 3);

  const auto extents_wout_halo = std::array { static_cast<std::size_t>(n_cells[0]),
                                              static_cast<std::size_t>(n_cells[1]),
                                              static_cast<std::size_t>(n_cells[2]) };
  const auto extents_with_halo =
    std::array { extents_wout_halo[0] + 2uz * emf::halo_size,
                 extents_wout_halo[1] + 2uz * emf::halo_size,
                 extents_wout_halo[2] + 2uz * emf::halo_size };

  auto construct_comm_buffs =
    [extents_wout_halo, extents_with_halo, sim](const auto id) {
      using EB_t = emf_comm_buffs::hollow_grid_EB;
      using J_t  = emf_comm_buffs::hollow_grid_J;
      sim.get().tiles.template emplace<emf_comm_buffs>(
        id,
        EB_t(extents_wout_halo),
        EB_t(extents_wout_halo),
        J_t(extents_with_halo));
    };

  return te::just() | te::then([sim, construct_comm_buffs] {
           auto where_to_construct = std::vector<simulation_context::tile_id_type> {};

           for(const auto& [id, _]:
               sim.get()
                 .tiles.view<virtual_tile_tag>(entt::exclude<emf_comm_buffs>)
                 .each()) {
             where_to_construct.push_back(id);
           }
           for(const auto& [id, _]:
               sim.get()
                 .tiles.view<boundary_tile_tag>(entt::exclude<emf_comm_buffs>)
                 .each()) {
             where_to_construct.push_back(id);
           }

           std::ranges::for_each(where_to_construct, construct_comm_buffs);
         });
}

auto
  update_send_buff_B(const std::reference_wrapper<simulation_context> sim)
{
  namespace te = tyvi::exec;
  return te::just() | te::continues_on(te::thread_pool_scheduler {}) | te::then([sim] {
           const auto w = tyvi::mdgrid_work {};

           for(auto&& [_, comm_buffs, yee, __]: sim.get()
                                                  .view_tiles<
                                                    emf_comm_buffs,
                                                    const emf::YeeLattice,
                                                    const boundary_tile_tag>()) {
             comm_buffs.B.set_from_mds(w, yee.nonhalo_submds(yee.mds_B()));
           }

           w.wait();
         });
}

auto
  comm_external_B(const std::reference_wrapper<simulation_context> sim)
{
  namespace te = tyvi::exec;

  namespace pmpi   = pika::mpi::experimental;
  using index_type = runko::cartesian_index<3>;
  auto recvs =
    te::just() | te::let_value([sim] {
      auto senders = std::vector<te::unique_any_sender<>> {};
      for(auto&& [_, idx, virt, comm_buffs]:
          sim.get()
            .view_tiles<const index_type, const virtual_tile_tag, emf_comm_buffs>()) {

        const auto span = comm_buffs.B.span();

        // Pedantic check, in case we want to change the data type.
        static_assert(std::same_as<emf_comm_buffs::value_type, float>);

        const auto rank = virt.source_info.rank;
        const auto tag  = virt.source_info.tag;

        senders.push_back(
          te::just(
            reinterpret_cast<void*>(span.data()),
            runko::checked_cast<int>(span.size()),
            MPI_FLOAT,
            rank,
            tag,
            MPI_COMM_WORLD)

          | te::continues_on(te::thread_pool_scheduler {}) |
          pmpi::transform_mpi(MPI_Irecv));
      }
      return te::when_all_vector(std::move(senders));
    });

  auto sends =
    te::just() | te::let_value([sim] {
      auto senders = std::vector<te::unique_any_sender<>> {};
      for(auto&& [_, idx, boundary, comm_buffs]:
          sim.get()
            .view_tiles<const index_type, const boundary_tile_tag, emf_comm_buffs>()) {

        const auto span = comm_buffs.B.span();

        // Pedantic check, in case we want to change the data type.
        static_assert(std::same_as<emf_comm_buffs::value_type, float>);

        for(const auto [rank, tag]: boundary.dest_infos) {
          senders.push_back(
            te::just(
              reinterpret_cast<void*>(span.data()),
              runko::checked_cast<int>(span.size()),
              MPI_FLOAT,
              rank,
              tag,
              MPI_COMM_WORLD) |
            te::continues_on(te::thread_pool_scheduler {}) |
            pmpi::transform_mpi(MPI_Isend));
        }
      }
      return te::when_all_vector(std::move(senders));
    });

  return te::schedule(te::thread_pool_scheduler {}) |
         te::let_value([sim] { return ensure_constructed_emf_comm_buffs(sim); }) |
         te::let_value([sim] { return update_send_buff_B(sim); }) |
         te::let_value([s = std::move(sends), r = std::move(recvs)] {
           return te::when_all(std::move(s), std::move(r));
         });
}
auto
  comm_local_B(const std::reference_wrapper<simulation_context> sim)
{

  namespace te = tyvi::exec;

  auto f = [sim] {
    for(const auto dir: moore_neigh_dirs<3>()) {
      for(auto&& [_, cart_neighs, yee, idx]: sim.get()
                                               .view_tiles<
                                                 cartesian_neighbors<3>,
                                                 emf::YeeLattice,
                                                 cartesian_index<3>,
                                                 local_tile_tag>()) {
        if(const auto neigh_id_opt = cart_neighs.get(dir)) {
          const auto neigh_id = neigh_id_opt.value();
          const auto dir_arr  = dir.to_vec<int>().data;
          if(const auto p = sim.get().tiles.try_get<emf::YeeLattice>(neigh_id)) {
            yee.set_B_in_subregion(dir_arr, *p);
          } else if(const auto p = sim.get().tiles.try_get<emf_comm_buffs>(neigh_id)) {
            const auto w = tyvi::mdgrid_work {};
            yee.set_B_in_subregion(w, dir_arr, p->B);
            w.wait();
          } else {
            throw std::logic_error("comm_local(emf_B): invalid neighbor");
          }
        } else {
          const auto d = dir.to_vec<int>();
          throw std::runtime_error { std::format(
            "comm_local(emf_B): {} {} {} has missing neighbor in direction {} {} {}",
            idx[0],
            idx[1],
            idx[2],
            d[0],
            d[1],
            d[2]) };
        }
      }
    }
  };
  return te::just() | te::then(f);
}

}  // namespace runko
