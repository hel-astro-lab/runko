#pragma once

#include "core/emf2/tile.h"
#include "core/particles_common.h"
#include "core/pic2/particle.h"
#include "corgi/corgi.h"
#include "corgi/tile.h"
#include "tools/config_parser.h"

#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

namespace pic2 {

namespace mpi = mpi4cpp::mpi;

/*! \brief PiC v2 tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf2::Tile
 */

template<std::size_t D>
class Tile : virtual public emf2::Tile<D>, virtual public corgi::Tile<D> {

  static_assert(D == 3);

  std::map<runko::particle, ParticleContainer> particle_buffs_;


public:
  /// The type in which pos, vel and weights are stored in.
  using value_type = ParticleContainer::value_type;

  /// Construct Tile based on the given config.
  ///
  /// FIXME: document how the initializeion is done?
  ///
  /// In addition to emf2::Tile ctor requirements,
  /// the given config has to contain values for:
  ///
  /// FIXME: what about photons (qp, mp)?
  ///
  /// `qx`: charge of x particle species (x in {e, i, p})
  /// `mx`: mass-to-charge of x particle species (x in {e, i, p})
  /// `delgam`: temperature(?)
  /// `temperature_ratio`: T_i / T_e
  /// `sigma`: magnetization number (omega_ce/omega_pe)^2, including gamma for inertia
  /// `c_omp`: simulation skin depth simulation
  /// `ppc`: particles per cell per species
  ///
  /// FIXME: figure out meaning, implement and document all options below:
  /// `npasses`: number of current filter passes
  /// `n_test_prtcls`: number of particles TestParticleWriter writes to disc
  /// `l0`: Nx * NxMesh / max_mode / c_omp (from pic-trubulence)
  /// `g0`: maximum attainable particle energy (from pic-trubulence)
  /// `t0`: time steps in units of light-crossing (from pic-trubulence)
  /// `{e,b,j,p}_norm`: default normalizations(?) (from pic-trubulence)
  ///
  /// Parameters for packet init:
  /// `zeta`: perturbation amplitude
  /// `ell`: packet width
  /// `impact_param`
  /// `two_wave[_reversed]`
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_idx,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;

  std::array<std::vector<value_type>, 3> get_positions(runko::particle);
  std::array<std::vector<value_type>, 3> get_velocities(runko::particle);
  std::vector<value_type> get_weights(runko::particle);
};


}  // namespace pic2
