#pragma once

#include "core/communication_common.h"
#include "core/emf2/tile.h"
#include "core/particles_common.h"
#include "core/pic2/particle.h"
#include "corgi/corgi.h"
#include "corgi/tile.h"
#include "pybind11/numpy.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

namespace pic2 {

namespace mpi = mpi4cpp::mpi;

enum class ParticlePusher { boris };
enum class FieldInterpolator { linear_1st };
enum class CurrentDepositer { zigzag_1st, zigzag_1st_atomic };

struct ParticleStateBatch {
  using container_type = std::array<pybind11::array_t<double>, 3>;

  container_type pos;
  container_type vel;
};

/*! \brief PiC v2 tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf2::Tile
 */

template<std::size_t D>
class Tile : virtual public emf2::Tile<D>, virtual public corgi::Tile<D> {

  static_assert(D == 3);

  std::map<std::size_t, ParticleContainer> particle_buffs_;
  ParticlePusher particle_pusher_;
  FieldInterpolator field_interpolator_;
  CurrentDepositer current_depositer_;
  std::map<std::size_t, std::size_t> amount_of_particles_to_be_received_;

  using subregion_particle_buff =
    std::map<std::array<int, 3>, std::map<std::size_t, ParticleContainer>>;

  std::optional<subregion_particle_buff> particles_in_subregion_ {};
  void split_particles_to_subregions();

  // Required for global periodic boundary for tiles.
  std::array<ParticleContainer::value_type, 3> global_coordinate_mins_;
  std::array<ParticleContainer::value_type, 3> global_coordinate_maxs_;

public:
  /// The type in which pos and vel are stored in.
  using value_type = ParticleContainer::value_type;

  /// Construct Tile based on the given config.
  ///
  /// FIXME: document how the initialization is done?
  ///
  /// In addition to emf2::Tile ctor requirements,
  /// the given config has to contain values for:
  ///
  /// `qx`: charge of x:th particle species (x is natural number)
  /// `mx`: mass-to-charge of x:th particle species (x is natural number)
  /// `delgam`: temperature(?)
  /// `temperature_ratio`: T_i / T_e
  /// `sigma`: magnetization number (omega_ce/omega_pe)^2, including gamma for inertia
  /// `c_omp`: simulation skin depth resolution
  /// `ppc`: particles per cell per species
  /// `particle_pusher`: scheme to update particles velocities and positions
  /// `fields_interpolator`: scheme to interpolate E and B fields to particles
  ///
  /// Note that particle charges and masses qx and mx are read in order: 0, 1, ...
  /// If a i:th mass and charge are missing, the search is stopped.
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

  std::array<std::vector<value_type>, 3> get_positions(std::size_t);
  std::array<std::vector<value_type>, 3> get_velocities(std::size_t);


  using particle_generator =
    std::function<std::vector<runko::ParticleState>(double, double, double)>;

  /// Inject particles based on given generator.
  ///
  /// Generator is called for each cell coordinates.
  ///
  /// Particle type is assumed to be configured.
  void inject_to_each_cell(std::size_t particle_type, particle_generator);

  /// Inject given particles.
  ///
  /// Particle type is assumed to be configured.
  void inject(std::size_t particle_type, const std::vector<runko::ParticleState>&);

  using batch_array = pybind11::array_t<double>;
  using batch_particle_generator =
    std::function<ParticleStateBatch(batch_array, batch_array, batch_array)>;

  /// Inject particles based on given generator.
  ///
  /// Generator is called once with all cell coordinates.
  ///
  /// Particle type is assumed to be configured.
  void batch_inject_to_cells(std::size_t particle_type, batch_particle_generator);

  /// Push particles updating their velocities and positions.
  void push_particles();

  /// Deposit current from all particls.
  void deposit_current();

  std::vector<mpi4cpp::mpi::request>
    send_data(mpi4cpp::mpi::communicator& /*comm*/, int dest, int mode, int tag)
      override;

  std::vector<mpi4cpp::mpi::request>
    recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag)
      override;

  void pairwise_moore_communication_prelude(const int) override;
  void pairwise_moore_communication_postlude(const int) override;

  /// Get particles from haloregion of the other with comm_mode::pic_particle.
  ///
  /// Forward other communication modes to emf2::Tile.
  /// Assumes that the other tile is pic2::Tile or its descendant.
  void pairwise_moore_communication(
    const corgi::Tile<D>& /* other */,
    const std::array<int, D> dir_to_other,
    const int /* mode */
    ) override;
};


}  // namespace pic2
