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

void
  push_particles(runko::simulation_context& sim)
{

  using yee_value_type = emf::YeeLattice::value_type;
  const auto particle_pusher =
    sim.get_n_set_config<pic::ParticlePusher>([](auto&& conf) {
      const auto p = conf.template get_or_throw<std::string>("particle_pusher");
      if(p == "boris") {
        return pic::ParticlePusher::boris;
      } else if(p == "higuera_cary") {
        return pic::ParticlePusher::higuera_cary;
      } else if(p == "faraday") {
        return pic::ParticlePusher::faraday;
      } else {
        const auto msg = std::format("{} is not supported particle pusher.", p);
        throw std::runtime_error { msg };
      }
    });

  const auto field_interpolator =
    sim.get_n_set_config<pic::FieldInterpolator>([](auto&& conf) {
      const auto p = conf.template get_or_throw<std::string>("field_interpolator");
      if(p == "linear_1st") {
        return pic::FieldInterpolator::linear_1st;
      } else if(p == "linear_1st_unrolled") {
        return pic::FieldInterpolator::linear_1st_unrolled;
      } else {
        const auto msg = std::format("{} is not supported field_interpolator.", p);
        throw std::runtime_error { msg };
      }
    });

  const auto cfl = sim.config.template get_or_throw<double>("cfl");
  for(auto&& [_, yee, particles, idx, id_gen]: sim.template view_tiles<
                                               emf::YeeLattice,
                                               pic::particle_containers,
                                               const runko::cartesian_index<3>,
                                               pic::particle_id_generator,
                                               runko::local_tile_tag>()) {
    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);
    const auto origo_pos =
      std::array { static_cast<yee_value_type>(gc.mins()[0]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[1]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[2]) - emf::halo_size };

    auto push_impl = [&](const auto& interpolator) {
      for(auto& [_, pbuff]: particles) {
        switch(particle_pusher) {
          case ParticlePusher::boris:
            pbuff.push_particles_boris(cfl, interpolator);
            break;
          case ParticlePusher::higuera_cary:
            pbuff.push_particles_higuera_cary(cfl, interpolator);
            break;
          case ParticlePusher::faraday:
            pbuff.push_particles_faraday(cfl, interpolator);
            break;
          default:
            throw std::logic_error {
              "pic::Tile::push_particles: unkown particle pusher"
            };
        }
      }
    };

    switch(field_interpolator) {
      case FieldInterpolator::linear_1st:
        push_impl(yee.interpolate_EB_linear_1st(origo_pos));
        break;
      case FieldInterpolator::linear_1st_unrolled:
        push_impl(yee.interpolate_EB_linear_1st_unrolled(origo_pos));
        break;
      default:
        throw std::logic_error {
          "pic::Tile::push_particles: unkown field interpolator"
        };
    }
  }
}

void
  deposit_current(runko::simulation_context& sim)
{
  const auto current_depositer =
    sim.get_n_set_config<pic::CurrentDepositer>([](auto&& conf) {
      const auto p = conf.template get_or_throw<std::string>("current_depositer");
      if(p == "zigzag" or p == "zigzag_1st") {
        return pic::CurrentDepositer::zigzag_1st;
      } else if(p == "zigzag_1st_atomic") {
        return pic::CurrentDepositer::zigzag_1st_atomic;
      } else {
        const auto msg = std::format("{} is not supported current depositer.", p);
        throw std::runtime_error { msg };
      }
    });
  const auto cfl = sim.config.template get_or_throw<double>("cfl");

  for(auto&& [tile_id, yee, particles, idx, id_gen]: sim.template view_tiles<
                                                     emf::YeeLattice,
                                                     pic::particle_containers,
                                                     const runko::cartesian_index<3>,
                                                     pic::particle_id_generator,
                                                     runko::local_tile_tag>()) {
    yee.clear_current();

    using yee_value_type = emf::YeeLattice::value_type;
    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);
    const auto origo_pos =
      std::array { static_cast<yee_value_type>(gc.mins()[0]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[1]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[2]) - emf::halo_size };


    switch(current_depositer) {
      case CurrentDepositer::zigzag_1st:
        for(const auto& [_, pcontainer]: particles) {
          yee.deposit_current(pcontainer.current_zigzag_1st(origo_pos, cfl));
        }
        break;
      case CurrentDepositer::zigzag_1st_atomic: {
        struct J_cache {
          using type = runko::VecGrid<emf::YeeLattice::value_type>;
          type cache;
        };

        auto cache_ptr = sim.tiles.try_get<J_cache>(tile_id);
        if(not cache_ptr) {
          cache_ptr = &sim.tiles.emplace<J_cache>(
            tile_id,
            J_cache::type(yee.extents_with_halo()));
        }
        auto& generated_J = cache_ptr->cache;

        const auto genJmds = generated_J.mds();
        tyvi::mdgrid_work {}
          .for_each_index(
            genJmds,
            [=](const auto idx, const auto tidx) { genJmds[idx][tidx] = 0; })
          .wait();

        for(const auto& [_, pcontainer]: particles) {
          pcontainer.current_zigzag_1st(generated_J, origo_pos, cfl);
        }

        yee.deposit_current(generated_J);
        break;
      }
      default:
        throw std::logic_error { "pic::deposit_current: unkown current depositer" };
    }

    if(const auto corrJ_ptr = sim.tiles.try_get<pic::correction_J>(tile_id)) {
      if(corrJ_ptr->pending) {
        yee.deposit_current(corrJ_ptr->J);
        const auto corrJmds = corrJ_ptr->J.mds();
        tyvi::mdgrid_work {}
          .for_each_index(
            corrJmds,
            [=](const auto idx, const auto tidx) { corrJmds[idx][tidx] = 0; })
          .wait();
        corrJ_ptr->pending = false;
      }
    }
  }
}

void
  sort_particles(runko::simulation_context& sim)
{
  for(auto&& [_, yee, particles, idx]: sim.template view_tiles<
                                       emf::YeeLattice,
                                       pic::particle_containers,
                                       const runko::cartesian_index<3>,
                                       runko::local_tile_tag>()) {

    const auto m = yee.grid_mapping_with_halo();
    using M      = decltype(m);

    using F = pic::ParticleContainer::value_type;

    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);
    const auto origo_pos = std::array { static_cast<F>(gc.mins()[0]) - emf::halo_size,
                                        static_cast<F>(gc.mins()[1]) - emf::halo_size,
                                        static_cast<F>(gc.mins()[2]) - emf::halo_size };
    using Vec3F          = toolbox::Vec3<F>;

    auto score = [=](const F x, const F y, const F z) {
      const auto dx  = Vec3F(x, y, z) - Vec3F(origo_pos);
      const auto idx = dx.template as<typename M::index_type>();

      return m(idx[0], idx[1], idx[2]);
    };

    for(auto& [_, pbuff]: particles) { pbuff.sort(score); }
  }
}


void
  register_reflector_wall(
    runko::simulation_context& sim,
    const pic::reflector_wall& wall)
{
  sim
    .get_n_set_config<pic::reflectors>(
      [](auto&&) { return pic::reflectors { .walls {} }; })
    .walls.push_back(wall);
}

void
  reflect_particles(runko::simulation_context& sim)
{

  if(not sim.has_config<pic::reflectors>()) { return; }
  auto& reflectors = sim.get_config<pic::reflectors>();
  if(reflectors.walls.empty()) { return; }

  const auto e = toolbox::get_extent_list(
    sim.config,
    "n_cells_per_tile",
    3,
    toolbox::with_halo<emf::halo_size>);
  const auto cfl = sim.config.template get_or_throw<double>("cfl");

  for(auto&& [tile_id, particles, idx]: sim.template view_tiles<
                                        pic::particle_containers,
                                        const runko::cartesian_index<3>>()) {

    const auto gc = runko::global_coordinates(sim, idx.template as<double>().data);

    const auto mins            = gc.mins();
    const auto maxs            = gc.maxs();
    const auto wall_is_in_tile = [&](const pic::reflector_wall& w) {
      using vt = pic::reflector_wall::value_type;
      return w.walloc >= static_cast<vt>(mins[0]) - cfl &&
             w.walloc <= static_cast<vt>(maxs[0]);
    };
    if(std::ranges::none_of(reflectors.walls, wall_is_in_tile)) { continue; }

    auto corr_J_ptr = sim.tiles.try_get<pic::correction_J>(tile_id);
    if(not corr_J_ptr) {
      corr_J_ptr = &sim.tiles.emplace<pic::correction_J>(
        tile_id,
        true,
        pic::correction_J::type(std::array { e[0], e[1], e[2] }));
    }

    corr_J_ptr->pending = true;

    using yee_value_type = emf::YeeLattice::value_type;
    const auto origo_pos =
      std::array { static_cast<yee_value_type>(gc.mins()[0]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[1]) - emf::halo_size,
                   static_cast<yee_value_type>(gc.mins()[2]) - emf::halo_size };

    for(const auto& wall: reflectors.walls) {
      if(!wall_is_in_tile(wall)) { continue; }

      for(auto& [_, pbuff]: particles) {
        pbuff.reflect_at_wall(wall, corr_J_ptr->J, origo_pos, cfl);
      }
    }
  }
}

void
  advance_reflector_walls(runko::simulation_context& sim)
{

  if(not sim.has_config<pic::reflectors>()) { return; }
  auto& reflectors = sim.get_config<pic::reflectors>();
  if(reflectors.walls.empty()) { return; }

  const auto cfl = sim.config.template get_or_throw<double>("cfl");

  for(auto& wall: reflectors.walls) {
    wall.walloc += wall.betawall * pic::reflector_wall::value_type(cfl);
  }
}

}  // namespace pic
