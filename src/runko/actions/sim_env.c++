// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/actions/args.h"
#include "runko/actions/emf.h"
#include "runko/actions/env.h"
#include "runko/actions/pic.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/comm/external.h"
#include "runko/comm/local.h"
#include "runko/communication_common.h"
#include "tyvi/actions_ast.h"
#include "tyvi/actions_list.h"

#include <functional>
#include <variant>

namespace runko {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace rn = std::ranges;
namespace rv = std::views;
namespace py = pybind11;

tyvi::actions::sexpr
  build_sim_env(simulation_context& sim)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  // Due to hipcc compiler bug, env not be non-const.
  // It would be better to have it be non-const and be moved into
  // std::visit(ta::list_append, ...) but this workaround propably
  // is not a performance killer even if we have to do some extra copies.
  const auto env = ta::list(
    ta::cons(
      runko::symbol::set_cartesian_neighbors,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<simulation_context>>(args) |
               te::then(&runko::set_cartesian_neighbors<3>) |
               te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::set_cartesian_comm_infos,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<simulation_context>>(args) |
               te::then(&runko::set_cartesian_comm_infos<3>) |
               te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::ensure_constructed_yee_lattices,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::then(&emf::ensure_constructed_yee_lattices) |
               te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::set_EBJ,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 py::function,
                 py::function,
                 py::function>(args) |
               te::let_value(
                 [](const auto sim, const auto& Eh, const auto& Bh, const auto& Jh) {
                   return emf::set_EBJ(
                     sim,
                     Eh.template cast<emf::vector_field_function>(),
                     Bh.template cast<emf::vector_field_function>(),
                     Jh.template cast<emf::vector_field_function>());
                 });
      } }),
    ta::cons(
      runko::symbol::batch_set_EBJ,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 py::function,
                 py::function,
                 py::function,
                 py::function,
                 py::function,
                 py::function,
                 py::function,
                 py::function,
                 py::function>(args) |
               te::let_value([](
                               const auto sim,
                               const auto& Exh,
                               const auto& Eyh,
                               const auto& Ezh,
                               const auto& Bxh,
                               const auto& Byh,
                               const auto& Bzh,
                               const auto& Jxh,
                               const auto& Jyh,
                               const auto& Jzh) {
                 return emf::batch_set_EBJ(
                   sim,
                   Exh.template cast<emf::batch_vector_field_function>(),
                   Eyh.template cast<emf::batch_vector_field_function>(),
                   Ezh.template cast<emf::batch_vector_field_function>(),
                   Bxh.template cast<emf::batch_vector_field_function>(),
                   Byh.template cast<emf::batch_vector_field_function>(),
                   Bzh.template cast<emf::batch_vector_field_function>(),
                   Jxh.template cast<emf::batch_vector_field_function>(),
                   Jyh.template cast<emf::batch_vector_field_function>(),
                   Jzh.template cast<emf::batch_vector_field_function>());
               });
      } }),
    ta::cons(
      runko::symbol::add_current,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::let_value(&emf::add_current);
      } }),
    ta::cons(
      runko::symbol::register_antenna,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 emf::antenna_mode>(args) |
               te::then(&emf::register_antenna) | te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::deposit_antenna_current,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::let_value(&emf::deposit_antenna_current);
      } }),
    ta::cons(
      runko::symbol::push_e,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::let_value(&emf::push_e);
      } }),
    ta::cons(
      runko::symbol::push_half_b,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::let_value(&emf::push_half_b);
      } }),
    ta::cons(
      runko::symbol::filter_current,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::let_value(&emf::filter_current);
      } }),
    ta::cons(
      runko::symbol::register_edge_bc,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 emf::edge_bc>(args) |
               te::then(&emf::register_edge_bc) | te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::apply_edge_bc,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 emf::edge_bc,
                 runko::comm_mode>(args) |
               te::let_value(&emf::apply_edge_bc);
      } }),
    ta::cons(
      runko::symbol::apply_edge_bcs,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 runko::comm_mode>(args) |
               te::let_value(&emf::apply_edge_bcs);
      } }),
    ta::cons(
      runko::symbol::clear_bcs,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::then([](runko::simulation_context& sim) {
                 sim.template get_config<emf::boundary_conditions>() =
                   emf::boundary_conditions {};
                 return ta::null;
               });
      } }),
    ta::cons(
      runko::symbol::ensure_constructed_particle_containers,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::then(&pic::ensure_constructed_particle_containers) |
               te::then([] { return ta::null; });
      } }),
    ta::cons(
      runko::symbol::inject_to_each_cell,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 long,
                 py::function>(args) |
               te::then([](const auto sim, const auto ptype, const auto& pgen) {
                 pic::inject_to_each_cell(
                   sim,
                   runko::checked_cast<std::size_t>(ptype),
                   pgen.template cast<pic::particle_generator>());
                 return ta::null;
               });
      } }),
    ta::cons(
      runko::symbol::inject,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 long,
                 py::list>(args) |
               te::then([](const auto sim, const auto ptype, const auto& particles) {
                 pic::inject(
                   sim,
                   runko::checked_cast<std::size_t>(ptype),
                   particles
                     .template cast<std::vector<runko::ParticleState<double>>>());
                 return ta::null;
               });
      } }),
    ta::cons(
      runko::symbol::batch_inject_to_cells,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 long,
                 py::function>(args) |
               te::then([](const auto sim, const auto ptype, const auto& pgen) {
                 pic::batch_inject_to_cells(
                   sim,
                   runko::checked_cast<std::size_t>(ptype),
                   pgen.template cast<pic::batch_particle_generator>());
                 return ta::null;
               });
      } }),
    ta::cons(
      runko::symbol::batch_inject_in_x_stripe,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<runko::simulation_context>,
                 long,
                 py::function,
                 double,
                 double>(args) |
               te::then([](
                          const auto sim,
                          const auto ptype,
                          const auto& pgen,
                          const auto a,
                          const auto b) {
                 pic::batch_inject_in_x_stripe(
                   sim,
                   runko::checked_cast<std::size_t>(ptype),
                   pgen.template cast<pic::batch_particle_generator>(),
                   a,
                   b);
                 return ta::null;
               });
      } }),
    ta::cons(
      runko::symbol::push_particles,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(
                 args) |
               te::then(&pic::push_particles) | te::then([] { return ta::null; });
      } }),
    ta::cons(runko::symbol::current_context, std::ref(sim)),
    ta::cons(
      runko::symbol::comm_local,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<simulation_context>,
                 runko::comm_mode>(args) |
               te::let_value(&runko::comm_local);
      } }),
    ta::cons(
      runko::symbol::comm_external,
      ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
        return parse_atom_args<
                 std::reference_wrapper<simulation_context>,
                 runko::comm_mode>(args) |
               te::let_value(&runko::comm_external);
      } }));

  return std::visit(ta::list_append, env, runko::build_std_env());
}

}  // namespace runko
