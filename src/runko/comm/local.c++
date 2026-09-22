// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/comm/local.h"

#include "runko/comm/cartesian_grid.h"
#include "runko/comm/emf.h"
#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "tyvi/execution.h"

#include <array>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <utility>


namespace runko {

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

}  // namespace runko
