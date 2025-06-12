#include "core/pic2/tile.h"

#include "core/pic2/particle.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <ranges>

namespace pic2 {

template<std::size_t D>
Tile<D>::Tile(const toolbox::ConfigParser& conf) :
  corgi::Tile<D>(),
  emf2::Tile<D>(conf),
  particles_per_cell_ { conf.get<std::size_t>("ppc").value() }
{

  const auto number_of_cells = std::reduce(
    this->yee_lattice_extents_.begin(),
    this->yee_lattice_extents_.end(),
    1uz,
    std::multiplies<std::size_t> {});

  const auto number_of_particles = number_of_cells * particles_per_cell_;


  // electrons, ions, photons(?)
  const auto particle_container_args = std::array {
    ParticleContainerArgs { .N      = number_of_particles,
                            .charge = -conf.get<double>("qe").value(),
                            .mass   = std::fabs(conf.get<double>("me").value()) },
    ParticleContainerArgs { .N      = number_of_particles,
                            .charge = -conf.get<double>("qi").value(),
                            .mass   = std::fabs(conf.get<double>("mi").value()) },
    ParticleContainerArgs { .N      = number_of_particles,
                            .charge = -conf.get<double>("qp").value(),
                            .mass   = std::fabs(conf.get<double>("mp").value()) }
  };


  const auto Nspecies = conf.get<std::size_t>("Nspecies").value();
  for(const auto& arg: particle_container_args | std::views::take(Nspecies)) {
    particle_buffs_.emplace_back(arg);
  }
}

}  // namespace pic2

template class pic2::Tile<3>;
