// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "mpi.h"
#include "pika/mpi.hpp"
#include "runko/comm/cartesian_grid.h"
#include "runko/simulation_context.h"
#include "tyvi/actions.h"
#include "tyvi/execution.h"

#include <format>
#include <print>
#include <vector>


namespace runko {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
inline auto
  dump_emf(simulation_context& sim)
{
  auto open_file = [&sim] {
    MPI_File fh;
    const auto file_path =
      sim.config.get<std::string>("io_outdir").value() + "/emf.test";

    std::ignore = MPI_File_open(
      MPI_COMM_WORLD,
      file_path.c_str(),
      MPI_MODE_CREATE | MPI_MODE_WRONLY,
      MPI_INFO_NULL,
      &fh);
    return fh;
  };

  auto close_file = [](MPI_File fh) {
    std::ignore = MPI_File_close(&fh);
    return ta::null;
  };

  return te::schedule(te::thread_pool_scheduler {}) | te::then(open_file) |
         te::let_value([&sim](MPI_File fh) {
           auto senders       = std::vector<te::unique_any_sender<>> {};
           const auto n_tiles = toolbox::get_extent_list(sim.config, "n_tiles", 3);
           const auto n_cells =
             toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);
           const auto cells = static_cast<std::size_t>(
             (n_cells[0] + 6) * (n_cells[1] + 6) * (n_cells[2] + 6));
           const auto bytes_per_tile = ((4uz * cells * 9uz) / 4096uz) * 4096uz;


           using M = std::layout_left::mapping<std::dextents<std::size_t, 3>>;
           const auto m =
             M(std::dextents<std::ptrdiff_t, 3>(n_tiles[0], n_tiles[1], n_tiles[2]));

           auto tile_view = sim.view_tiles<
             const runko::cartesian_index<3>,
             const runko::local_tile_tag,
             const emf::YeeLattice>();

           const auto num_tiles =
             static_cast<std::size_t>(std::ranges::distance(tile_view));
           MPI_File_set_size(fh, static_cast<MPI_Offset>(num_tiles * bytes_per_tile));

           for(auto&& [_, idx, yee]: tile_view) {

             const auto N           = m(idx[0], idx[1], idx[2]);
             const auto tile_offset = N * bytes_per_tile;

             const auto E = yee.span_E();
             const auto B = yee.span_B();
             const auto J = yee.span_J();

             if(E.size_bytes() != 3uz * cells * sizeof(float)) {
               throw std::logic_error {
                 std::format("E {} {}", E.size_bytes(), 3uz * cells * sizeof(float))
               };
             }
             if(B.size_bytes() != 3uz * cells * sizeof(float)) {
               throw std::logic_error {
                 std::format("B {} {}", B.size_bytes(), 3uz * cells * sizeof(float))
               };
             }
             if(J.size_bytes() != 3uz * cells * sizeof(float)) {
               throw std::logic_error {
                 std::format("J {} {}", J.size_bytes(), 3uz * cells * sizeof(float))
               };
             }

             namespace pmpi = pika::mpi::experimental;
             senders.push_back(
               te::just(
                 fh,
                 tile_offset,
                 reinterpret_cast<const void*>(E.data()),
                 E.size(),
                 MPI_FLOAT) |
               pmpi::transform_mpi(MPI_File_iwrite_at));
             senders.push_back(
               te::just(
                 fh,
                 tile_offset + 3uz * cells * sizeof(float),
                 reinterpret_cast<const void*>(B.data()),
                 B.size(),
                 MPI_FLOAT) |
               pmpi::transform_mpi(MPI_File_iwrite_at));
             senders.push_back(
               te::just(
                 fh,
                 tile_offset + 2uz * 3uz * cells * sizeof(float),
                 reinterpret_cast<const void*>(J.data()),
                 J.size(),
                 MPI_FLOAT) |
               pmpi::transform_mpi(MPI_File_iwrite_at));
           }

           return te::when_all_vector(std::move(senders)) |
                  te::then([fh] { return fh; });
         }) |
         te::then(close_file);
}

}  // namespace runko
