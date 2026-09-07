// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/emf.h"
#include "runko/actions/env.h"
#include "runko/comm/emf.h"
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

tyvi::actions::sexpr
  build_sim_env(simulation_context& sim)
{
  namespace ta = tyvi::actions;
  namespace te = tyvi::exec;

  auto comm_local = [](const ta::sexpr& args) -> ta::sexpr_sender {
    const auto [sim, mode] = args_to_sim_n_comm_mode(args);
    return te::just() | te::continues_on(te::thread_pool_scheduler {}) |
           te::let_value([sim, mode]() -> ta::sexpr_sender {
             switch(mode) {
               case comm_mode::emf_B:
                 return comm_local_B(sim) | te::then([] { return ta::null; });
               default:
                 throw std::runtime_error { "comm_local: unregonized comm mode" };
             }

             return te::just(ta::null);
           });
  };

  auto comm_external = [](const ta::sexpr& args) -> ta::sexpr_sender {
    const auto [sim, mode] = args_to_sim_n_comm_mode(args);
    return te::just() | te::continues_on(te::thread_pool_scheduler {}) |
           te::let_value([sim, mode]() -> ta::sexpr_sender {
             switch(mode) {
               case comm_mode::emf_B:
                 return comm_external_B(sim) | te::then([] { return ta::null; });
               default:
                 throw std::runtime_error { "comm_external: unregonized comm mode" };
             }

             return te::just(ta::null);
           });
  };

  // Due to hipcc compiler bug, env not be non-const.
  // It would be better to have it be non-const and be moved into
  // std::visit(ta::list_append, ...) but this workaround propably
  // is not a performance killer even if we have to do some extra copies.
  const auto env = ta::list(
    ta::cons(runko::symbol::set_EBJ, ta::procedure { &runko::set_EBJ }),
    ta::cons(runko::symbol::current_context, std::ref(sim)),
    ta::cons(runko::symbol::comm_local, ta::procedure { comm_local }),
    ta::cons(runko::symbol::comm_external, ta::procedure { comm_external }));

  return std::visit(ta::list_append, env, runko::build_std_env());
}

}  // namespace runko
