// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/comm/emf.h"

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/comm/cartesian_grid.h"
#include "runko/comm/external.h"
#include "runko/communication_common.h"
#include "runko/coords.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "tyvi/execution.h"

#include <array>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <utility>

namespace emf {

void
  ensure_constructed_emf_comm_buffs(runko::simulation_context& sim)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  const auto n_cells = toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);

  const auto extents_wout_halo = std::array { static_cast<std::size_t>(n_cells[0]),
                                              static_cast<std::size_t>(n_cells[1]),
                                              static_cast<std::size_t>(n_cells[2]) };
  const auto extents_with_halo =
    std::array { extents_wout_halo[0] + 2uz * emf::halo_size,
                 extents_wout_halo[1] + 2uz * emf::halo_size,
                 extents_wout_halo[2] + 2uz * emf::halo_size };

  auto construct_comm_buffs =
    [extents_wout_halo, extents_with_halo, &sim](const auto id) {
      using EB_t = emf::comm_buffs::hollow_grid_EB;
      using J_t  = emf::comm_buffs::hollow_grid_J;
      sim.tiles.template emplace<emf::comm_buffs>(
        id,
        EB_t(extents_wout_halo),
        EB_t(extents_wout_halo),
        J_t(extents_with_halo));
    };

  auto where_to_construct = std::vector<runko::simulation_context::tile_id_type> {};

  for(const auto& [id, _]:
      sim.tiles.view<runko::virtual_tile_tag>(entt::exclude<emf::comm_buffs>).each()) {
    where_to_construct.push_back(id);
  }
  for(const auto& [id, _]:
      sim.tiles.view<runko::boundary_tile_tag>(entt::exclude<emf::comm_buffs>).each()) {
    where_to_construct.push_back(id);
  }

  std::ranges::for_each(where_to_construct, construct_comm_buffs);
}

void
  update_emf_send_buff(runko::simulation_context& sim, const runko::comm_mode mode)
{
  const auto w = tyvi::mdgrid_work {};

  for(auto&& [_, comm_buffs, yee, boundary]: sim.view_tiles<
                                             emf::comm_buffs,
                                             const emf::YeeLattice,
                                             const runko::boundary_tile_tag>()) {
    std::ignore = boundary;

    switch(mode) {
      case runko::comm_mode::emf_E:

        for(auto&& [_, comm_buffs, yee, __]: sim.view_tiles<
                                             emf::comm_buffs,
                                             const emf::YeeLattice,
                                             const runko::boundary_tile_tag>()) {
          comm_buffs.E.set_from_mds(w, yee.nonhalo_submds(yee.mds_E()));
        }
        w.wait();
        break;
      case runko::comm_mode::emf_B:

        for(auto&& [_, comm_buffs, yee, __]: sim.view_tiles<
                                             emf::comm_buffs,
                                             const emf::YeeLattice,
                                             const runko::boundary_tile_tag>()) {
          comm_buffs.B.set_from_mds(w, yee.nonhalo_submds(yee.mds_B()));
        }
        w.wait();
        break;
      case runko::comm_mode::emf_J:

        for(auto&& [_, comm_buffs, yee, __]: sim.view_tiles<
                                             emf::comm_buffs,
                                             const emf::YeeLattice,
                                             const runko::boundary_tile_tag>()) {
          comm_buffs.J.set_from_mds(w, yee.mds_J());
        }
        w.wait();
        break;
      default:
        throw std::logic_error {
          std::format("update_emf_send_buff({}): unregonized comm_mode", mode)
        };
    }
  }
}


tyvi::actions::sexpr_sender
  comm_external(runko::simulation_context& x, const runko::comm_mode mode)
{
  namespace te = tyvi::exec;

  namespace pmpi   = pika::mpi::experimental;
  using index_type = runko::cartesian_index<3>;


  auto get_span = [](emf::comm_buffs& buffs, const runko::comm_mode mode) {
    switch(mode) {
      case runko::comm_mode::emf_E: return buffs.E.span(); break;
      case runko::comm_mode::emf_B: return buffs.B.span(); break;
      case runko::comm_mode::emf_J: return buffs.J.span(); break;
      default:
        throw std::runtime_error {
          std::format("emf::comm_external: unregonized comm_mode: {}", mode)
        };
    }
  };

  const auto comm = [mode, sim = std::ref(x)] {
    switch(mode) {
      case runko::comm_mode::emf_E:
        return runko::make_communicator<runko::comm_mode::emf_E>(sim).comm;
        break;
      case runko::comm_mode::emf_B:
        return runko::make_communicator<runko::comm_mode::emf_B>(sim).comm;
        break;
      case runko::comm_mode::emf_J:
        return runko::make_communicator<runko::comm_mode::emf_J>(sim).comm;
        break;
      default:
        throw std::runtime_error {
          std::format("emf::comm_external: unregonized comm_mode: {}", mode)
        };
    }
  }();

  auto recvs = [mode, get_span, comm](runko::simulation_context& sim) {
    auto senders = std::vector<te::unique_any_sender<>> {};
    for(auto&& [_, idx, virt, comm_buffs]: sim.view_tiles<
                                           const index_type,
                                           const runko::virtual_tile_tag,
                                           emf::comm_buffs>()) {

      const auto span = get_span(comm_buffs, mode);

      // Pedantic check, in case we want to change the data type.
      static_assert(std::same_as<emf::comm_buffs::value_type, float>);

      const auto rank = virt.source_info.rank;
      const auto tag  = virt.source_info.tag;

      senders.push_back(
        te::just(
          reinterpret_cast<void*>(span.data()),
          runko::checked_cast<int>(span.size()),
          MPI_FLOAT,
          rank,
          tag,
          comm)

        | te::continues_on(te::thread_pool_scheduler {}) |
        pmpi::transform_mpi(MPI_Irecv));
    }
    return te::when_all_vector(std::move(senders));
  };

  auto sends = [mode, get_span, comm](runko::simulation_context& sim) {
    auto senders = std::vector<te::unique_any_sender<>> {};
    for(auto&& [_, idx, boundary, comm_buffs]: sim.view_tiles<
                                               const index_type,
                                               const runko::boundary_tile_tag,
                                               emf::comm_buffs>()) {

      const auto span = get_span(comm_buffs, mode);

      // Pedantic check, in case we want to change the data type.
      static_assert(std::same_as<emf::comm_buffs::value_type, float>);

      for(const auto [rank, tag]: boundary.dest_infos) {
        senders.push_back(
          te::just(
            reinterpret_cast<void*>(span.data()),
            runko::checked_cast<int>(span.size()),
            MPI_FLOAT,
            rank,
            tag,
            comm) |
          te::continues_on(te::thread_pool_scheduler {}) |
          pmpi::transform_mpi(MPI_Isend));
      }
    }
    return te::when_all_vector(std::move(senders));
  };

  return te::just() |
         te::then(std::bind_front(&ensure_constructed_emf_comm_buffs, std::ref(x))) |
         te::then(std::bind_front(&update_emf_send_buff, std::ref(x), mode)) |
         te::let_value([sim = std::ref(x), sends, recvs] {
           return te::when_all(sends(sim), recvs(sim));
         }) |
         te::drop_value() | te::then([] { return tyvi::actions::null; });
}


tyvi::actions::sexpr_sender
  comm_local(runko::simulation_context& sim, const runko::comm_mode mode)
{
  namespace te = tyvi::exec;

  auto f = [sim = std::ref(sim), mode] {
    for(const auto dir: runko::moore_neigh_dirs<3>()) {

      const auto w = tyvi::mdgrid_work {};
      for(auto&& [_, cart_neighs, yee, idx]: sim.get()
                                               .view_tiles<
                                                 runko::cartesian_neighbors<3>,
                                                 emf::YeeLattice,
                                                 runko::cartesian_index<3>,
                                                 runko::local_tile_tag>()) {

        const auto neigh_id_opt = cart_neighs.get(dir);
        if(not neigh_id_opt) {
          const auto d = dir.to_vec<int>();
          throw std::runtime_error { std::format(
            "comm_local({}): {} {} {} has missing neighbor in direction {} {} {}",
            mode,
            idx[0],
            idx[1],
            idx[2],
            d[0],
            d[1],
            d[2]) };
        }

        const auto neigh_id = neigh_id_opt.value();
        const auto dir_arr  = dir.to_vec<int>().data;

        if(const auto p = sim.get().tiles.try_get<emf::YeeLattice>(neigh_id)) {
          switch(mode) {
            case runko::comm_mode::emf_E: yee.set_E_in_subregion(w, dir_arr, *p); break;
            case runko::comm_mode::emf_B: yee.set_B_in_subregion(w, dir_arr, *p); break;
            case runko::comm_mode::emf_J: yee.set_J_in_subregion(w, dir_arr, *p); break;
            case runko::comm_mode::emf_J_exchange:
              yee.add_to_J_from_subregion(w, dir_arr, *p);
              break;
            default:
              throw std::logic_error { std::format(
                "comm_local({}): unhandled comm_mode for non-virtual neighbor",
                mode) };
          }
        } else if(const auto p = sim.get().tiles.try_get<emf::comm_buffs>(neigh_id)) {

          switch(mode) {
            case runko::comm_mode::emf_E:
              yee.set_E_in_subregion(w, dir_arr, p->E);
              break;
            case runko::comm_mode::emf_B:
              yee.set_B_in_subregion(w, dir_arr, p->B);
              break;
            case runko::comm_mode::emf_J:
              yee.set_J_in_subregion(w, dir_arr, p->J);
              break;
            case runko::comm_mode::emf_J_exchange:
              yee.add_to_J_from_subregion(w, dir_arr, p->J);
              break;
            default:
              throw std::logic_error { std::format(
                "comm_local({}): unhandled comm_mode for virtual neighbor",
                mode) };
          }
        } else {
          throw std::logic_error(std::format("comm_local({}): invalid neighbor", mode));
        }
      }

      w.wait();
    }

    return tyvi::actions::null;
  };
  return te::just() | te::then(f);
}

}  // namespace emf
