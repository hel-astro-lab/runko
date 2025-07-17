#include "core/pic2/tile.h"

#include "core/particles_common.h"
#include "core/pic2/particle.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <ranges>
#include <string>

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
      // Negation charges and abs(m) are as they are in pytools/pic/*
      return pic2::ParticleContainerArgs { .N = 0, .charge = -q, .mass = std::fabs(m) };
    };

    return conf.get<double>(q_label)
      .and_then(q_to_qm_tuple)
      .transform(qm_tuple_to_args);
  };

  if(const auto args = make_opt_args("qe", "me")) {
    where_to.insert_or_assign(
      runko::particle::electron,
      pic2::ParticleContainer { args.value() });
  }

  if(const auto args = make_opt_args("qi", "mi")) {
    where_to.insert_or_assign(
      runko::particle::ion,
      pic2::ParticleContainer { args.value() });
  }

  if(const auto args = make_opt_args("qp", "mp")) {
    where_to.insert_or_assign(
      runko::particle::photon,
      pic2::ParticleContainer { args.value() });
  }
}

}  // namespace

namespace pic2 {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_grid_idx,
  const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf2::Tile<D>(tile_grid_idx, conf)
{
  construct_particle_buffs(particle_buffs_, conf);
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_positions(const runko::particle p)
{
  return particle_buffs_.at(p).get_positions();
}

template<std::size_t D>
std::array<std::vector<typename Tile<D>::value_type>, 3>
  Tile<D>::get_velocities(const runko::particle p)
{
  return particle_buffs_.at(p).get_velocities();
}

template<std::size_t D>
std::vector<typename Tile<D>::value_type>
  Tile<D>::get_weights(const runko::particle p)
{
  return particle_buffs_.at(p).get_weights();
}

template<std::size_t D>
void
  Tile<D>::inject_to_each_cell(
    const runko::particle particle_type,
    particle_generator pgen)
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
    for(const auto p: pgen(x, y, z)) {
        new_particles.push_back(p);
    }
  }

  particle_buffs_.at(particle_type).add_particles(new_particles);
}

}  // namespace pic2

template class pic2::Tile<3>;
