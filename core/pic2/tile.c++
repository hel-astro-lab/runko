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

}  // namespace pic2

template class pic2::Tile<3>;
