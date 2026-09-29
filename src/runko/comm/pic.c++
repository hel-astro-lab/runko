// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/comm/pic.h"

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/actions/pic.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/comm/external.h"
#include "runko/coords.h"
#include "thrust/memory.h"

#include <ranges>

namespace pic {

void
  ensure_constructed_pic_comm_buffs(runko::simulation_context& sim)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  pic::ensure_constructed_particle_containers(sim);

  for(const auto& [id, _]:
      sim.tiles.view<runko::virtual_tile_tag>(entt::exclude<pic::external_comm_buffs>)
        .each()) {
    sim.tiles.template emplace<pic::external_comm_buffs>(id);
  }
  for(const auto& [id]:
      sim.tiles.view<runko::local_tile_tag>(entt::exclude<pic::external_comm_buffs>)
        .each()) {
    sim.tiles.template emplace<pic::external_comm_buffs>(id);
  }

  for(const auto& [id]:
      sim.tiles.view<runko::local_tile_tag>(entt::exclude<pic::local_comm_buffs>)
        .each()) {
    std::ignore = sim.tiles.template emplace<pic::local_comm_buffs>(id);
  }
}

void
  pack_outgoing_particles_for_tile(
    const runko::global_coordinates_closure<3> gc,
    pic::external_comm_buffs& comm_buffs,
    pic::particle_containers& pcontainers)
{

  const auto mins  = runko::cast_array<pic::ParticleContainer::value_type>(gc.mins());
  const auto maxs  = runko::cast_array<pic::ParticleContainer::value_type>(gc.maxs());
  const auto x_div = std::array { mins[0], maxs[0] };
  const auto y_div = std::array { mins[1], maxs[1] };
  const auto z_div = std::array { mins[2], maxs[2] };

  comm_buffs.subregion_particle_ends.resize(27 * pcontainers.size());
  comm_buffs.subregion_particle_buff.resize(0);

  // ptypes are assumed always to be contiguous 0, 1, ..., N
  for(auto& [ptype, pcontainer]: pcontainers) {
    auto spans = pcontainer.divide_to_subregions(
      comm_buffs.subregion_particle_buff,
      x_div,
      y_div,
      z_div);

    for(const auto& [dir, span]: spans) {
      comm_buffs.subregion_particle_ends.at(27 * ptype + dir.neighbor_index()) =
        std::get<1>(span);
    }
  }
}

auto
  pack_outgoing_particles(runko::simulation_context& sim)
{
  namespace te = tyvi::exec;
  auto senders = std::vector<te::unique_any_sender<>> {};
  for(auto&& [_, idx, particles, comm_buffs]: sim.view_tiles<
                                              const runko::cartesian_index<3>,
                                              pic::particle_containers,
                                              pic::external_comm_buffs,
                                              runko::local_tile_tag>()) {
    senders.push_back(
      te::just(
        runko::global_coordinates(sim, idx.template as<double>().data),
        std::ref(comm_buffs),
        std::ref(particles)) |
      te::continues_on(te::thread_pool_scheduler {}) |
      te::then(&pack_outgoing_particles_for_tile));
  }
  return te::when_all_vector(std::move(senders));
}

tyvi::actions::sexpr_sender
  comm_external_particles(runko::simulation_context& x)
{
  namespace te     = tyvi::exec;
  namespace pmpi   = pika::mpi::experimental;
  using index_type = runko::cartesian_index<3>;

  const auto number_of_ptypes = [&] {
    for(auto i = 0uz; true; ++i) {
      const auto q_label = std::format("q{}", i);
      const auto m_label = std::format("m{}", i);

      if(x.config.contains(q_label) and x.config.contains(m_label)) {
        continue;
      } else {
        return i;
      }
    }
  }();


  const auto comm = runko::make_communicator<runko::comm_mode::pic_particle>(x).comm;
  using value_type =
    decltype(pic::external_comm_buffs::subregion_particle_buff)::value_type;

  auto virt_tiles = [comm, number_of_ptypes](runko::simulation_context& sim) {
    auto senders = std::vector<te::unique_any_sender<>> {};
    for(auto&& [id, idx, virt, comm_buffs]: sim.view_tiles<
                                            const index_type,
                                            const runko::virtual_tile_tag,
                                            pic::external_comm_buffs>()) {

      const auto rank = virt.source_info.rank;
      const auto tag  = virt.source_info.tag;

      senders.push_back(
        te::schedule(te::thread_pool_scheduler {}) |
        te::let_value([rank, tag, comm, number_of_ptypes, &comm_buffs, id] {
          comm_buffs.subregion_particle_ends.resize(number_of_ptypes * 27);
          std::ignore = id;
          return te::just(
                   reinterpret_cast<void*>(comm_buffs.subregion_particle_ends.data()),
                   runko::checked_cast<int>(comm_buffs.subregion_particle_ends.size()),
                   MPI_UINT64_T,
                   rank,
                   tag,
                   comm) |
                 pmpi::transform_mpi(MPI_Irecv);
        }) |
        te::let_value([rank, tag, comm, &comm_buffs] {
          comm_buffs.subregion_particle_buff.resize(
            comm_buffs.subregion_particle_ends.back());
          // This is a bit pendantic :)
          static_assert(CHAR_BIT == 8uz);
          return te::just(
                   reinterpret_cast<void*>(thrust::raw_pointer_cast(
                     comm_buffs.subregion_particle_buff.data())),
                   runko::checked_cast<int>(
                     comm_buffs.subregion_particle_buff.size() * sizeof(value_type)),
                   MPI_BYTE,
                   rank,
                   tag,
                   comm) |
                 pmpi::transform_mpi(MPI_Irecv);
        }));
    }
    return te::when_all_vector(std::move(senders));
  };

  auto boundary_tiles = [comm](runko::simulation_context& sim) {
    auto senders = std::vector<te::unique_any_sender<>> {};
    for(auto&& [_, idx, boundary, comm_buffs]: sim.view_tiles<
                                               const index_type,
                                               const runko::boundary_tile_tag,
                                               pic::external_comm_buffs>()) {

      for(const auto [rank, tag]: boundary.dest_infos) {
        senders.push_back(
          te::schedule(te::thread_pool_scheduler {}) |
          te::let_value([rank, tag, comm, &comm_buffs] {
            return te::just(
                     reinterpret_cast<void*>(comm_buffs.subregion_particle_ends.data()),
                     runko::checked_cast<int>(
                       comm_buffs.subregion_particle_ends.size()),
                     MPI_UINT64_T,
                     rank,
                     tag,
                     comm) |
                   pmpi::transform_mpi(MPI_Isend);
          }) |
          te::let_value([rank, tag, comm, &comm_buffs] {
            // This is a bit pendantic :)
            static_assert(CHAR_BIT == 8uz);
            return te::just(
                     reinterpret_cast<void*>(thrust::raw_pointer_cast(
                       comm_buffs.subregion_particle_buff.data())),
                     runko::checked_cast<int>(
                       comm_buffs.subregion_particle_buff.size() * sizeof(value_type)),
                     MPI_BYTE,
                     rank,
                     tag,
                     comm) |
                   pmpi::transform_mpi(MPI_Isend);
          }));
      }
    }
    return te::when_all_vector(std::move(senders));
  };

  return te::just() |
         te::then(std::bind_front(&ensure_constructed_pic_comm_buffs, std::ref(x))) |
         te::let_value(std::bind_front(&pack_outgoing_particles, std::ref(x))) |
         te::let_value([sim = std::ref(x), virt_tiles, boundary_tiles] {
           return te::when_all(virt_tiles(sim), boundary_tiles(sim));
         }) |
         te::then([] { return tyvi::actions::null; });
}

tyvi::actions::sexpr_sender
  comm_local_particles(runko::simulation_context& x)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  const auto number_of_ptypes = [&] {
    for(auto i = 0uz; true; ++i) {
      const auto q_label = std::format("q{}", i);
      const auto m_label = std::format("m{}", i);

      if(x.config.contains(q_label) and x.config.contains(m_label)) {
        continue;
      } else {
        return i;
      }
    }
  }();

  auto prep_local_comm_buffs = [number_of_ptypes](
                                 const runko::simulation_context& sim,
                                 const runko::cartesian_index<3>& idx,
                                 const runko::cartesian_neighbors<3>& cart_neighs,
                                 pic::local_comm_buffs& buffs) {
    for(const auto dir: runko::moore_neigh_dirs<3>()) {
      const auto neigh_id_opt = cart_neighs.get(dir);

      if(not neigh_id_opt) {
        const auto d = dir.to_vec<int>();
        throw std::runtime_error { std::format(
          "comm_local_particles(): {} {} {} has missing neighbor in direction {} {} "
          "{}",
          idx[0],
          idx[1],
          idx[2],
          d[0],
          d[1],
          d[2]) };
      }

      const auto neigh_id  = neigh_id_opt.value();
      const auto neigh_ptr = sim.tiles.try_get<pic::external_comm_buffs>(neigh_id);

      if(not neigh_ptr) {
        throw std::logic_error("comm_local_particle(): invalid neighbor");
      }

      const auto inverted_dir = dir.inverted();

      namespace rv = std::views;
      for(const auto ptype: rv::iota(0uz, number_of_ptypes)) {
        const auto& ends = neigh_ptr->subregion_particle_ends;
        const auto& buff = neigh_ptr->subregion_particle_buff;

        const auto index = 27 * ptype + inverted_dir.neighbor_index();
        const auto end   = static_cast<std::ptrdiff_t>(ends.at(index));
        const auto begin =
          static_cast<std::ptrdiff_t>(index == 0 ? 0uz : ends.at(index - 1));
        const auto p = thrust::raw_pointer_cast(buff.data());

        buffs[ptype].push_back(
          std::span<const runko::ParticleState<pic::ParticleContainer::value_type>>(
            std::ranges::next(p, begin),
            std::ranges::next(p, end)));
      }
    }
  };

  const auto [gce_mins, gce_maxs] = runko::global_coordinate_extents<3>(x);

  auto gather_from_others =
    [gce_mins,
     gce_maxs](pic::local_comm_buffs& buffs, pic::particle_containers& pcontainers) {
      for(const auto& [ptype, spans]: buffs) {
        using T = pic::ParticleContainer::value_type;
        pcontainers.at(ptype).append(
          spans,
          runko::cast_array<T>(gce_mins),
          runko::cast_array<T>(gce_maxs));
      }
      buffs.clear();
    };

  auto f = [prep_local_comm_buffs, gather_from_others](runko::simulation_context& sim) {
    auto senders = std::vector<te::unique_any_sender<>> {};

    for(auto&& [_, cart_neighs, particles, idx, comm_buffs]:
        sim.view_tiles<
          runko::cartesian_neighbors<3>,
          pic::particle_containers,
          runko::cartesian_index<3>,
          pic::local_comm_buffs,
          runko::local_tile_tag>()) {
      senders.push_back(
        te::just(std::cref(sim), idx, cart_neighs, std::ref(comm_buffs)) |
        te::continues_on(te::thread_pool_scheduler {}) |
        te::then(prep_local_comm_buffs) |
        te::then(
          std::bind_front(
            gather_from_others,
            std::ref(comm_buffs),
            std::ref(particles))));
    }

    return te::when_all_vector(std::move(senders));
  };
  return te::just(std::ref(x)) | te::let_value(f) | te::then([] { return ta::null; });
}

}  // namespace pic
