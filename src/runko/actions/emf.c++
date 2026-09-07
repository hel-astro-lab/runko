// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/emf.h"

#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "runko/coords.h"
#include "runko/emf/common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/config_parser.h"
#include "tyvi/execution.h"

#include <algorithm>
#include <functional>
#include <print>
#include <tuple>

namespace runko {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace py = pybind11;

void
  ensure_constructed_yee_lattices(const std::reference_wrapper<simulation_context> sim)
{
  auto where_to_construct = std::vector<runko::simulation_context::tile_id_type> {};

  for(auto&& id:
      sim.get().tiles.view<runko::local_tile_tag>(entt::exclude<emf::YeeLattice>)) {
    where_to_construct.push_back(id);
  }
  if(where_to_construct.empty()) { return; }

  const auto n_cells =
    toolbox::get_extent_list(sim.get().config, "n_cells_per_tile", 3);
  const auto args = emf::make_yee_ctor_args_from_vector(n_cells);
  std::ranges::for_each(where_to_construct, [&](const auto x) {
    sim.get().tiles.template emplace<emf::YeeLattice>(x, args);
  });
}

ta::sexpr_sender
  set_EBJ(const ta::sexpr& args)
{
  const auto args0 = std::get<ta::cons>(args);
  const auto arg0  = std::get<ta::atom>(args0.car());
  const auto args1 = std::get<ta::cons>(args0.cdr());
  const auto arg1  = std::get<ta::atom>(args1.car());
  const auto args2 = std::get<ta::cons>(args1.cdr());
  const auto arg2  = std::get<ta::atom>(args2.car());
  const auto args3 = std::get<ta::cons>(args2.cdr());
  const auto arg3  = std::get<ta::atom>(args3.car());

  const auto sim_opt =
    ta::atom_cast<std::reference_wrapper<runko::simulation_context>>(arg0);
  if(not sim_opt) {
    throw std::runtime_error("set_EBJ not invoked with simulation_context");
  }
  const auto sim = sim_opt.value();

  auto Eopt = ta::atom_cast<py::function>(arg1);
  auto Bopt = ta::atom_cast<py::function>(arg2);
  auto Jopt = ta::atom_cast<py::function>(arg3);

  static constexpr auto err_msg =
    "set_EBJ not invoked with three functions: (float, float, float) -> (float, "
    "float, float)";

  if(not(Eopt and Bopt and Jopt)) {
    std::println("foo");
    throw std::runtime_error { err_msg };
  }

  const auto Eh = Eopt.value();
  const auto Bh = Bopt.value();
  const auto Jh = Jopt.value();

  using F = std::function<std::tuple<double, double, double>(double, double, double)>;

  auto E = Eh.cast<F>();
  auto B = Bh.cast<F>();
  auto J = Jh.cast<F>();

  return te::just(sim) | te::then(&ensure_constructed_yee_lattices) |
         te::then(
           [sim, E = std::move(E), B = std::move(B), J = std::move(J)]() -> ta::sexpr {
             for(auto&& [id, yee, idx]: sim.get()
                                          .view_tiles<
                                            emf::YeeLattice,
                                            const runko::cartesian_index<3>,
                                            runko::local_tile_tag>()) {
               const auto gc =
                 runko::global_coordinates(sim.get(), idx.template as<double>().data);
               auto f =
                 [&](const std::size_t i, const std::size_t j, const std::size_t k) {
                   const auto [x, y, z] = gc(i, j, k);
                   const auto ex        = E(x + 0.5, y, z);
                   const auto ey        = E(x, y + 0.5, z);
                   const auto ez        = E(x, y, z + 0.5);
                   const auto bx        = B(x, y + 0.5, z + 0.5);
                   const auto by        = B(x + 0.5, y, z + 0.5);
                   const auto bz        = B(x + 0.5, y + 0.5, z);
                   const auto jx        = J(x + 0.5, y, z);
                   const auto jy        = J(x, y + 0.5, z);
                   const auto jz        = J(x, y, z + 0.5);

                   return emf::YeeLatticeFieldsAtPoint { .Ex = std::get<0>(ex),
                                                         .Ey = std::get<1>(ey),
                                                         .Ez = std::get<2>(ez),
                                                         .Bx = std::get<0>(bx),
                                                         .By = std::get<1>(by),
                                                         .Bz = std::get<2>(bz),
                                                         .Jx = std::get<0>(jx),
                                                         .Jy = std::get<1>(jy),
                                                         .Jz = std::get<2>(jz) };
                 };
               yee.set_EBJ(f);
             }
             return ta::null;
           });
}


}  // namespace runko
