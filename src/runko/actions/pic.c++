// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/pic.h"

#include "runko/comm/cartesian_grid.h"
#include "runko/coords.h"
#include "runko/pic/particle.h"
#include "runko/tools/math.h"
#include "tyvi/execution.h"
#include "tyvi/mdspan.h"

#include <algorithm>
#include <array>
#include <optional>
#include <string>
#include <vector>

namespace pic {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;

class [[nodiscard]] particle_id_generator {
  std::vector<std::size_t> ordinals_;
  runko::prtc_id_type tag_;

public:
  explicit constexpr particle_id_generator(const runko::prtc_id_type tag) : tag_ { tag }
  {
    if(tag >= 2uz << 24uz) {
      throw std::runtime_error {
        "particle_id_generator requires tag to be less than 2^24."
      };
    }
  }

  [[nodiscard]] constexpr runko::prtc_id_type consume_next_id(const std::size_t ptype)
  {
    if(this->ordinals_.size() <= ptype) { ordinals_.resize(ptype + 1uz); }
    const auto ordinal = static_cast<runko::prtc_id_type>(this->ordinals_.at(ptype)++);

    if(ordinal >= 2uz << 40uz) {
      throw std::runtime_error {
        "PIC tile ran out of particle ids (2^40 per species)!"
      };
    }

    return (this->tag_ << 40uz) | ordinal;
  }
};

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


  const auto n_tiles = toolbox::get_extent_list(sim.config, "n_tiles", 3);


  auto handle_tile = [&](const auto x) {
    // First construct particle_id_generator.

    const auto idx_p = sim.tiles.template try_get<runko::cartesian_index<3>>(x);
    if(idx_p == nullptr) {
      throw std::runtime_error {
        "Initializing pic data to tile requires runko::cartesian_index<3>."
      };
    }
    const auto idx = *idx_p;
    const auto tag = std::layout_right::mapping {
      std::dextents<std::size_t, 3> { n_tiles[0], n_tiles[1], n_tiles[2] }
    }(idx[0], idx[1], idx[2]);

    std::ignore = sim.tiles.template emplace_or_replace<pic::particle_id_generator>(
      x,
      static_cast<runko::prtc_id_type>(tag));

    // Now construcct the particle containers.

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


void
  inject_to_each_cell(
    runko::simulation_context& sim,
    const std::size_t particle_type,
    pic::particle_generator pgen)
{
  const auto e = toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);

  for(auto&& [_, particles, idx, id_gen]: sim.template view_tiles<
                                          pic::particle_containers,
                                          const runko::cartesian_index<3>,
                                          pic::particle_id_generator,
                                          runko::local_tile_tag>()) {

    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);
    std::vector<runko::ParticleState<double>> new_particles {};

    for(const auto [i, j, k]:
        tyvi::sstd::index_space(std::mdspan((int*)nullptr, e[0], e[1], e[2]))) {
      const auto [x, y, z] = gc(i, j, k);
      for(const auto p: pgen(x, y, z)) {
        new_particles.push_back(p);
        new_particles.back().id = id_gen.consume_next_id(particle_type);
      }
    }

    particles.at(particle_type).add_particles(new_particles);
  }
}

void
  inject(
    runko::simulation_context& sim,
    std::size_t particle_type,
    std::vector<runko::ParticleState<double>> new_particles)
{
  for(auto&& [_, particles, idx, id_gen]: sim.template view_tiles<
                                          pic::particle_containers,
                                          const runko::cartesian_index<3>,
                                          pic::particle_id_generator,
                                          runko::local_tile_tag>()) {
    for(auto& p: new_particles) { p.id = id_gen.consume_next_id(particle_type); }
    particles.at(particle_type).add_particles(new_particles);
  }
}

void
  batch_inject_to_cells(
    runko::simulation_context& sim,
    std::size_t particle_type,
    batch_particle_generator pgen)
{
  const auto gce = runko::global_coordinate_extents<3>(sim);
  pic::batch_inject_in_x_stripe(
    sim,
    particle_type,
    std::move(pgen),
    gce[0][0],
    gce[1][0]);
}

void
  batch_inject_in_x_stripe(
    runko::simulation_context& sim,
    std::size_t particle_type,
    batch_particle_generator pgen,
    const double x_left,
    const double x_right)
{

  const auto e = toolbox::get_extent_list(sim.config, "n_cells_per_tile", 3);

  for(auto&& [_, particles, idx, id_gen]: sim.template view_tiles<
                                          pic::particle_containers,
                                          const runko::cartesian_index<3>,
                                          pic::particle_id_generator,
                                          runko::local_tile_tag>()) {
    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);

    const auto tile_xmin = static_cast<double>(gc.mins()[0]);
    const auto tile_xmax = static_cast<double>(gc.maxs()[0]);

    // early exit if stripe does not overlap this tile
    if(x_right <= tile_xmin || x_left >= tile_xmax) continue;


    // find x-cells overlapping the stripe (O(1) since dx=1)
    // cell i spans [tile_xmin+i, tile_xmin+i+1)
    const auto diff_left  = x_left - tile_xmin;
    const auto diff_right = x_right - tile_xmin;

    const auto i_begin = diff_left <= 0.0
                           ? std::size_t { 0 }
                           : std::ranges::min(
                               static_cast<std::size_t>(e[0]),
                               static_cast<std::size_t>(sstd::floor(diff_left)));

    const auto i_end = diff_right >= static_cast<double>(e[0])
                         ? static_cast<std::size_t>(e[0])
                         : static_cast<std::size_t>(sstd::ceil(diff_right));

    if(i_begin >= i_end) continue;

    const auto nx_stripe = i_end - i_begin;
    const auto cells_in_total =
      std::array { nx_stripe * static_cast<std::size_t>(e[1] * e[2]) };
    auto x = pybind11::array_t<double>(cells_in_total);
    auto y = pybind11::array_t<double>(cells_in_total);
    auto z = pybind11::array_t<double>(cells_in_total);

    {
      auto xv = x.template mutable_unchecked<1>();
      auto yv = y.template mutable_unchecked<1>();
      auto zv = z.template mutable_unchecked<1>();

      std::size_t n = 0;
      for(std::size_t i = i_begin; i < i_end; ++i) {
        for(std::size_t j = 0; j < static_cast<std::size_t>(e[1]); ++j) {
          for(std::size_t k = 0; k < static_cast<std::size_t>(e[2]); ++k) {
            const auto c = gc(i, j, k);
            xv(n)        = c[0];
            yv(n)        = c[1];
            zv(n)        = c[2];
            ++n;
          }
        }
      }
    }

    const auto state_batch = pgen(x, y, z);
    const auto batch_size  = static_cast<std::size_t>(state_batch.pos[0].shape(0));

    {
      const auto assert_shapes = [&](const auto& a) {
        if(a.ndim() != 1) {
          throw std::runtime_error {
            "pic::Tile::batch_inject_in_x_stripe: given batch must be one dimensional."
          };
        }

        if(static_cast<std::size_t>(a.shape(0)) != batch_size) {
          throw std::runtime_error {
            "pic::Tile::batch_inject_in_x_stripe: batches must have same length."
          };
        }
      };

      assert_shapes(state_batch.pos[0]);
      assert_shapes(state_batch.pos[1]);
      assert_shapes(state_batch.pos[2]);
      assert_shapes(state_batch.vel[0]);
      assert_shapes(state_batch.vel[1]);
      assert_shapes(state_batch.vel[2]);
    }

    const auto posx_view = state_batch.pos[0].template unchecked<1>();
    const auto posy_view = state_batch.pos[1].template unchecked<1>();
    const auto posz_view = state_batch.pos[2].template unchecked<1>();

    const auto velx_view = state_batch.vel[0].template unchecked<1>();
    const auto vely_view = state_batch.vel[1].template unchecked<1>();
    const auto velz_view = state_batch.vel[2].template unchecked<1>();

    auto states = std::vector<runko::ParticleState<double>> {};
    states.reserve(batch_size);
    for(const auto n: std::views::iota(0uz, batch_size)) {
      // partial edge cells are generated whole; keep only the part inside the stripe
      if(posx_view(n) < x_left or posx_view(n) >= x_right) { continue; }
      states.push_back(
        runko::ParticleState<double> {
          .pos = { posx_view(n), posy_view(n), posz_view(n) },
          .vel = { velx_view(n), vely_view(n), velz_view(n) },
          .id  = id_gen.consume_next_id(particle_type) });
    }

    particles.at(particle_type).add_particles(states);
  }
}

}  // namespace pic
