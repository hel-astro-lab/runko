#pragma once

#include "core/emf/tile.h"
#include "corgi/tile.h"
#include "mpi.h"

#include <array>
#include <concepts>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>

namespace {


/// Used to unify implementation.
///
/// Assumes that f returns a double.
void
  write_average_field_value(
    std::invocable<const emf::Tile<3>&> auto&& f,
    const runko::size_t lap,
    const std::string& path,
    corgi::Grid<3>& grid)
{
  namespace rn = std::ranges;
  namespace rv = std::ranges::views;

  const auto local_tile_indices = grid.get_local_tiles();

  const auto num_local_tiles = rn::size(local_tile_indices);
  if(num_local_tiles == 0u) {
    throw std::runtime_error {
      "write_average_field_value assumes that every rank has at least one tile."
    };
  }


  const auto emf_tiles =
    local_tile_indices | rv::transform([&](const auto cid) -> emf::Tile<3>& {
      try {
        return dynamic_cast<emf::Tile<3>&>(grid.get_tile(cid));
      } catch(const std::bad_cast& ex) {
        throw std::runtime_error { std::format(
          "write_average_field_value assumes that all tiles are"
          "emf::Tiles or its descendants. Orginal exception: {}",
          ex.what()) };
      }
    });


  auto total_value = std::invoke([&]() -> double {
    const auto local_values = emf_tiles | rv::transform(std::forward<decltype(f)>(f));
    return std::reduce(rn::begin(local_values), rn::end(local_values));
  });

  auto total_area = std::invoke([&]() -> double {
    const auto local_areas = emf_tiles | rv::transform([](const auto& tile) {
                               const auto [nx, ny, nz] = tile.extents_wout_halo();
                               return nx * ny * nz;
                             });

    return static_cast<double>(
      std::reduce(rn::begin(local_areas), rn::end(local_areas)));
  });


  std::array<MPI_Request, 2> requests;

  if(grid.comm.rank() == 0) {
    MPI_Ireduce(
      MPI_IN_PLACE,
      &total_value,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      MPI_COMM_WORLD,
      &requests[0]);
  } else {
    MPI_Ireduce(
      &total_value,
      nullptr,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      MPI_COMM_WORLD,
      &requests[0]);
  }

  if(grid.comm.rank() == 0) {
    MPI_Ireduce(
      MPI_IN_PLACE,
      &total_area,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      MPI_COMM_WORLD,
      &requests[1]);
  } else {
    MPI_Ireduce(
      &total_area,
      nullptr,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      MPI_COMM_WORLD,
      &requests[1]);
  }

  std::array<MPI_Status, 2> statuses;
  if(MPI_SUCCESS != MPI_Waitall(requests.size(), requests.data(), statuses.data())) {
    throw std::runtime_error { "MPI reduction in write_average_field_value failed!" };
  }

  if(grid.comm.rank() == 0) {
    std::fstream file { path, std::ios::out | std::ios::app };
    file << lap << ' ' << total_value / total_area << '\n';
  }
}
}  // namespace


namespace emf {

void
  write_average_B_energy_density(
    const runko::size_t lap,
    const std::string path,
    corgi::Grid<3>& grid)
{
  write_average_field_value(
    [](const auto& tile) { return tile.total_energy_B(); },
    lap,
    path,
    grid);
}

void
  write_average_E_energy_density(
    const runko::size_t lap,
    const std::string path,
    corgi::Grid<3>& grid)
{
  write_average_field_value(
    [](const auto& tile) { return tile.total_energy_E(); },
    lap,
    path,
    grid);
}

}  // namespace emf
