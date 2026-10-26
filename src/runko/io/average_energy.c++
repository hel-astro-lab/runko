// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/io/average_energy.h"

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/actions/pic.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/io/fs_utils.h"
#include "runko/io/snapshot.h"
#include "runko/particles_common.h"
#include "thrust/count.h"
#include "thrust/iterator/counting_iterator.h"
#include "tyvi/execution.h"
#include "tyvi/mdgrid.h"

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace runko {

namespace ta = tyvi::actions;
namespace te = tyvi::exec;

namespace {

int
  comm_world_rank(const std::source_location location = std::source_location::current())
{
  int comm_rank = -1;
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank)) {
    throw std::runtime_error {
      std::format("{}:{}: MPI_Comm_rank failed", location.file_name(), location.line())
    };
  }
  return comm_rank;
}

template<typename T>
auto
  reduce_sum_to_root(
    const std::span<T> buff,
    const MPI_Datatype type,
    const int comm_rank)
{
  const auto n   = checked_cast<int>(buff.size());
  namespace pmpi = pika::mpi::experimental;
  return te::just(
           comm_rank == 0 ? MPI_IN_PLACE : buff.data(),
           buff.data(),
           n,
           type,
           MPI_SUM,
           0,
           MPI_COMM_WORLD) |
         pmpi::transform_mpi(MPI_Ireduce);
}

/// Number of alive particles in the container.
std::size_t
  count_alive(const pic::ParticleContainer& container, const tyvi::mdgrid_work& w)
{
  const auto slots_begin = thrust::counting_iterator<runko::index_t>(0uz);

  return static_cast<std::size_t>(thrust::count_if(
    w.on_this(),
    slots_begin,
    slots_begin + container.ssize(),
    [ids_mds = container.ids_mds()](const runko::index_t n) {
      return ids_mds[n][] != runko::dead_prtc_id;
    }));
}

/// Append line to {outdir}/{filename}.
auto
  append_line(
    const std::string& outdir,
    const std::string_view filename,
    const std::string_view line)
{
  std::filesystem::create_directories(outdir);
  auto file = std::ofstream { std::format("{}/{}", outdir, filename), std::ios::app };
  file << line << '\n';
}

/// Used to unify implementation of field energy densities.
tyvi::actions::sexpr_sender
  average_field_energy_density(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir_arg,
    const std::string_view filename,
    std::invocable<const emf::YeeLattice&> auto f)
{
  return te::just(std::array<double, 2> {}) |
         te::continues_on(te::thread_pool_scheduler {}) |
         te::let_value(
           [&sim, lap, outdir_arg = std::move(outdir_arg), filename, f](auto& data) {
             double& energy = data[0];
             double& volume = data[1];

             for(auto&& [_, yee]:
                 sim.view_tiles<const emf::YeeLattice, runko::local_tile_tag>()) {
               const auto [nx, ny, nz] = yee.extents_wout_halo();
               energy += std::invoke(f, yee);
               volume += static_cast<double>(nx * ny * nz);
             }

             const auto comm_rank = comm_world_rank();

             return reduce_sum_to_root(std::span<double>(data), MPI_DOUBLE, comm_rank) |
                    te::then([&energy,
                              &volume,
                              comm_rank,
                              &sim,
                              outdir_arg = std::move(outdir_arg),
                              filename,
                              lap] {
                      if(comm_rank == 0) {
                        const auto avg = volume == 0.0 ? 0.0 : energy / volume;
                        append_line(
                          resolve_outdir(sim.config, outdir_arg),
                          filename,
                          std::format("{} {}", lap, avg));
                      }

                      return ta::null;
                    });
           });
}
}  // namespace

tyvi::actions::sexpr_sender
  average_B_energy_density(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir)
{
  return average_field_energy_density(
    sim,
    lap,
    std::move(outdir),
    "average_B_energy_density.txt",
    &emf::YeeLattice::total_energy_B);
}

tyvi::actions::sexpr_sender
  average_E_energy_density(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir)
{
  return average_field_energy_density(
    sim,
    lap,
    std::move(outdir),
    "average_E_energy_density.txt",
    &emf::YeeLattice::total_energy_E);
}

tyvi::actions::sexpr_sender
  average_kinetic_energy(
    runko::simulation_context& sim,
    const long lap,
    std::optional<std::string> outdir_arg)
{
  static constexpr auto filename = std::string_view { "average_kinetic_energy.txt" };

  // Number of species has to be same on every rank for the reductions,
  // so it is taken from the config instead of the particle containers.
  const auto nspecies = count_species(sim.config);

  return te::just(
           std::vector<double>(nspecies, 0.0),
           std::vector<std::uint64_t>(nspecies, 0)) |
         te::continues_on(te::thread_pool_scheduler {}) |
         te::let_value([&sim, lap, outdir_arg = std::move(outdir_arg), nspecies](
                         auto& total_energies,
                         auto& total_counts) {
           const auto w = tyvi::mdgrid_work {};

           for(auto&& [_, particles]:
               sim
                 .view_tiles<const pic::particle_containers, runko::local_tile_tag>()) {
             for(const auto& [species, container]: particles) {
               if(species >= nspecies) {
                 throw std::runtime_error { std::format(
                   "average_kinetic_energy: particle species {} is not configured.",
                   species) };
               }
               total_energies[species] += container.total_kinetic_energy();
               total_counts[species] += count_alive(container, w);
             }
           }

           const auto comm_rank = comm_world_rank();


           // Make sure that all ranks execute the reduces in order.
           return reduce_sum_to_root(
                    std::span<double>(total_energies),
                    MPI_DOUBLE,
                    comm_rank) |
                  te::let_value(
                    [s = std::span<std::uint64_t>(total_counts), comm_rank] {
                      return reduce_sum_to_root(s, MPI_UINT64_T, comm_rank);
                    }) |
                  te::then([comm_rank,
                            nspecies,
                            &total_energies,
                            &total_counts,
                            &sim,
                            lap,
                            outdir_arg = std::move(outdir_arg)] {
                    if(comm_rank == 0) {
                      auto line = std::format("{}", lap);
                      for(auto s = 0uz; s < nspecies; ++s) {
                        const auto avg =
                          total_counts[s] == 0
                            ? 0.0
                            : total_energies[s] / static_cast<double>(total_counts[s]);
                        line += std::format(" {}", avg);
                      }

                      append_line(
                        resolve_outdir(sim.config, outdir_arg),
                        filename,
                        line);
                    }
                    return ta::null;
                  });
         });
}

}  // namespace runko
