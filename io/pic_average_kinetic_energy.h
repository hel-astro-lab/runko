#pragma once

#include "core/pic/tile.h"
#include "corgi/tile.h"
#include "mpi.h"

#include <concepts>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace pic {

void
  write_average_kinetic_energy(
    const std::size_t lap,
    const std::string path,
    corgi::Grid<3>& grid)
{
  namespace rn = std::ranges;
  namespace rv = std::ranges::views;

  const auto local_tile_indices = grid.get_local_tiles();

  const auto num_local_tiles = std::ranges::size(local_tile_indices);
  if(num_local_tiles == 0u) {
    throw std::runtime_error {
      "write_average_kinetic_energy assumes that every rank has at least one tile."
    };
  }


  auto pic_tiles =
    local_tile_indices | rv::transform([&](const auto cid) -> pic::Tile<3>& {
      try {
        return dynamic_cast<pic::Tile<3>&>(grid.get_tile(cid));
      } catch(const std::bad_cast& ex) {
        throw std::runtime_error { std::format(
          "average_kinetic_energy assumes that all tiles are"
          "pic::Tiles or its descendants. Orginal exception: {}",
          ex.what()) };
      }
    });

  auto local_total_energy_n_number_of_particles = [&](const std::size_t particle_type) {
    auto total_energies = pic_tiles | rv::transform([&](auto& tile) -> double {
                            return tile.total_kinetic_energy(particle_type);
                          });

    auto number_of_particles = pic_tiles | rv::transform([&](auto& tile) -> double {
                                 return tile.number_of_particles(particle_type);
                               });

    return std::tuple {
      std::reduce(total_energies.begin(), total_energies.end()),
      std::reduce(number_of_particles.begin(), number_of_particles.end())
    };
  };

  const auto num_species   = (*rn::begin(pic_tiles)).number_of_species();
  auto total_energies      = std::vector<double>(num_species);
  auto number_of_particles = std::vector<std::size_t>(num_species);

  auto requests = std::vector<MPI_Request> {};

  for(const auto n: rv::iota(0u, num_species)) {
    const auto [E, N]      = local_total_energy_n_number_of_particles(n);
    total_energies[n]      = E;
    number_of_particles[n] = N;

    requests.emplace_back();

    if(grid.comm.rank() == 0) {
      MPI_Ireduce(
        MPI_IN_PLACE,
        &total_energies[n],
        1,
        MPI_DOUBLE,
        MPI_SUM,
        0,
        MPI_COMM_WORLD,
        &requests.back());
    } else {
      MPI_Ireduce(
        &total_energies[n],
        nullptr,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        0,
        MPI_COMM_WORLD,
        &requests.back());
    }

    static_assert(sizeof(std::size_t) == sizeof(std::uint64_t));

    requests.emplace_back();

    if(grid.comm.rank() == 0) {
      MPI_Ireduce(
        MPI_IN_PLACE,
        &number_of_particles[n],
        1,
        MPI_UINT64_T,
        MPI_SUM,
        0,
        MPI_COMM_WORLD,
        &requests.back());
    } else {
      MPI_Ireduce(
        &number_of_particles[n],
        nullptr,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        0,
        MPI_COMM_WORLD,
        &requests.back());
    }
  }

  auto statuses = std::vector<MPI_Status>(requests.size());
  if(MPI_SUCCESS != MPI_Waitall(requests.size(), requests.data(), statuses.data())) {
    throw std::runtime_error {
      "MPI reduction while calculating average energies failed!"
    };
  }

  if(grid.comm.rank() == 0) {
    std::fstream file { path, std::ios::out | std::ios::app };
    file << lap << ' ';
    for(const auto n: rv::iota(0u, num_species)) {
      const auto e =
        number_of_particles[n] == 0u
          ? 0.0
          : total_energies[n] / static_cast<double>(number_of_particles[n]);
      file << e << ' ';
    }
    file << "\n";
  }
}
}  // namespace pic
