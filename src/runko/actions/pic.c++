// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/pic.h"

#include "runko/comm/cartesian_grid.h"
#include "runko/pic/particle.h"
#include "tyvi/execution.h"
#include "tyvi/mdspan.h"

#include <optional>
#include <string>

namespace pic {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;

void
  ensure_constructed_particle_containers(runko::simulation_context& sim)
{
  auto where_to_construct = std::vector<runko::simulation_context::tile_id_type> {};

  for(auto&& id:
      sim.tiles.view<runko::local_tile_tag>(entt::exclude<pic::particle_containers>)) {
    where_to_construct.push_back(id);
  }
  if(where_to_construct.empty()) { return; }

  using opt_args_t = std::optional<pic::ParticleContainerArgs>;

  // Return optional containing the arguments if q and m are found in conf.
  const auto make_opt_args =
    [&](const std::string& q_label, const std::string& m_label) -> opt_args_t {
    const auto q_to_qm_tuple = [&](const auto q) {
      return sim.config.get<double>(m_label).transform(
        [&](const auto m) { return std::tuple { q, m }; });
    };

    const auto qm_tuple_to_args = [&](const auto qm) {
      const auto [q, m] = qm;

      return pic::ParticleContainerArgs { .N = 0, .charge = q, .mass = m };
    };

    return sim.config.get<double>(q_label)
      .and_then(q_to_qm_tuple)
      .transform(qm_tuple_to_args);
  };

  const auto prealloc_per_species =
    sim.config.get<std::size_t>("prealloc_per_species").value_or(0uz);

  auto handle_tile = [&](const auto x) {
    auto& where_to = sim.tiles.template emplace<pic::particle_containers>(x);

    for(auto i = 0uz; true; ++i) {
      const auto q_label = std::format("q{}", i);
      const auto m_label = std::format("m{}", i);

      const auto pcontainer_args = make_opt_args(q_label, m_label);
      if(not pcontainer_args) { break; }

      auto container = pic::ParticleContainer { pcontainer_args.value() };
      if(prealloc_per_species > 0uz) { container.prealloc_dead(prealloc_per_species); }

      // insert_or_assign and not operator[], because if element is missing,
      // then operator[] will default construct it.
      // pic::ParticleContainer is not default constructible.
      std::ignore = where_to.insert_or_assign(i, std::move(container));
    }
  };

  std::ranges::for_each(where_to_construct, handle_tile);
}

}  // namespace pic
