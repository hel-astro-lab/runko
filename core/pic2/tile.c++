#include "core/pic2/tile.h"

#include "core/emf2/yee_lattice.h"
#include "core/particles_common.h"
#include "core/pic2/particle.h"

#include <cmath>
#include <format>
#include <functional>
#include <numeric>
#include <ranges>
#include <string>
#include <string_view>

namespace {

void
  construct_particle_buffs(auto& where_to, const toolbox::ConfigParser& conf)
{
  using opt_args_t = std::optional<pic2::ParticleContainerArgs>;

  // Return optional containing the arguments if q and m are found in conf.
  const auto make_opt_args =
    [&](const std::string& q_label, const std::string& m_label) -> opt_args_t {
    const auto q_to_qm_tuple = [&](const auto q) {
      return conf.get<double>(m_label).transform(
        [&](const auto m) { return std::tuple { q, m }; });
    };

    const auto qm_tuple_to_args = [&](const auto qm) {
      const auto [q, m] = qm;

      // Due to effects of density function to particle number,
      // we do not yet know how many particles there are going to be.
      return pic2::ParticleContainerArgs { .N = 0, .charge = q, .mass = m };
    };

    return conf.get<double>(q_label)
      .and_then(q_to_qm_tuple)
      .transform(qm_tuple_to_args);
  };

  for(auto i = 0uz; true; ++i) {
    const auto q_label = std::format("q{}", i);
    const auto m_label = std::format("m{}", i);

    if(const auto pcontainer_args = make_opt_args(q_label, m_label)) {
      // insert_or_assign and not operator[], because if element is missing,
      // then operator[] will default construct it.
      // pic2::ParticleContainer is not default constructible.
      std::ignore = where_to.insert_or_assign(
        i,
        pic2::ParticleContainer { pcontainer_args.value() });

    } else {
      break;
    }
  }
}

pic2::ParticlePusher
  parse_particle_pusher(const std::string_view p)
{
  if(p == "boris") {
    return pic2::ParticlePusher::boris;
  } else {
    const auto msg = std::format("{} is not supported particle pusher.", p);
    throw std::runtime_error { msg };
  }
}

pic2::FieldInterpolator
  parse_field_interpolator(const std::string_view p)
{
  if(p == "linear_1st") {
    return pic2::FieldInterpolator::linear_1st;
  } else {
    const auto msg = std::format("{} is not supported field_interpolator.", p);
    throw std::runtime_error { msg };
  }
}

pic2::CurrentDepositer
  parse_current_depositer(const std::string_view p)
{
  if(p == "zigzag" or p == "zigzag_1st") {
    return pic2::CurrentDepositer::zigzag_1st;
  } else {
    const auto msg = std::format("{} is not supported current depositer.", p);
    throw std::runtime_error { msg };
  }
}

}  // namespace

namespace pic2 {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_grid_idx,
  const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf2::Tile<D>(tile_grid_idx, conf),
  particle_pusher_ { parse_particle_pusher(
    conf.get_or_throw<std::string>("particle_pusher")) },
  field_interpolator_ { parse_field_interpolator(
    conf.get_or_throw<std::string>("field_interpolator")) },
  current_depositer_ { parse_current_depositer(
    conf.get_or_throw<std::string>("current_depositer")) }
{
  construct_particle_buffs(particle_buffs_, conf);
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_positions(const std::size_t p)
{
  return particle_buffs_.at(p).get_positions();
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_velocities(const std::size_t p)
{
  return particle_buffs_.at(p).get_velocities();
}

template<std::size_t D>
void
  Tile<D>::inject_to_each_cell(const std::size_t particle_type, particle_generator pgen)
{
  if(this->mins == this->maxs) {
    throw std::logic_error {
      "failed pic2::Tile::inject_to_each_cell precondition: corgi::tile::{mins,maxs} "
      "aren't set"
    };
  }

  /// FIXME: unify global coordinates from here and emf2::Tile::set_EBJ to mixin.
  const auto Lx = static_cast<double>(this->maxs[0] - this->mins[0]);
  const auto Ly = static_cast<double>(this->maxs[1] - this->mins[1]);
  const auto Lz = static_cast<double>(this->maxs[2] - this->mins[2]);

  const auto e = this->extents_wout_halo();

  auto global_coordinates = [&](const auto i, const auto j, const auto k) {
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(e[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(e[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(e[2]);

    return std::tuple<double, double, double> {
      static_cast<double>(this->mins[0]) + x_coeff * Lx,
      static_cast<double>(this->mins[1]) + y_coeff * Ly,
      static_cast<double>(this->mins[2]) + z_coeff * Lz
    };
  };

  std::vector<runko::ParticleState> new_particles {};

  namespace rv     = std::views;
  auto index_space = rv::cartesian_product(
    rv::iota(0uz, e[0]),
    rv::iota(0uz, e[1]),
    rv::iota(0uz, e[2]));
  for(const auto [i, j, k]: index_space) {
    const auto [x, y, z] = global_coordinates(i, j, k);
    for(const auto p: pgen(x, y, z)) { new_particles.push_back(p); }
  }

  particle_buffs_.at(particle_type).add_particles(new_particles);
}

template<std::size_t D>
void
  Tile<D>::push_particles(const std::size_t p)
{
  using yee_value_type = emf2::YeeLattice::value_type;
  const auto origo_pos =
    std::array { static_cast<yee_value_type>(this->mins[0]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[1]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[2]) - this->halo_size };

  pic2::ParticleContainer::InterpolatedEB_function ipol_func {};

  switch(field_interpolator_) {
    case FieldInterpolator::linear_1st: {
      ipol_func = std::bind_front(
        &emf2::YeeLattice::interpolate_EB_linear_1st,
        std::cref(this->yee_lattice_),
        origo_pos);
      break;
    }
    default:
      throw std::logic_error { "pic2::Tile::push_particles: unkown interpolator" };
  }


  switch(particle_pusher_) {
    case ParticlePusher::boris:
      particle_buffs_.at(p).push_particles_boris(this->cfl_, std::move(ipol_func));
      break;
    default:
      throw std::logic_error { "pic2::Tile::push_particles: unkown particle pusher" };
  }
}

template<std::size_t D>
void
  Tile<D>::deposit_current()
{
  this->yee_lattice_.clear_current();

  using yee_value_type = emf2::YeeLattice::value_type;
  const auto origo_pos =
    std::array { static_cast<yee_value_type>(this->mins[0]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[1]) - this->halo_size,
                 static_cast<yee_value_type>(this->mins[2]) - this->halo_size };

  switch(current_depositer_) {
    case CurrentDepositer::zigzag_1st:
      for(const auto& [_, pcontainer]: this->particle_buffs_) {
        this->yee_lattice_.deposit_current(
          pcontainer.current_zigzag_1st(origo_pos, this->cfl_));
      }
      break;
    default:
      throw std::logic_error { "pic2::Tile::push_particles: unkown current depositer" };
  }
}

}  // namespace pic2

template class pic2::Tile<3>;
