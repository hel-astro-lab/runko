#pragma once

#include "core/communication_common.h"
#include "core/emf/tile.h"
#include "core/particles_common.h"
#include "core/pic/particle.h"
#include "core/pic/reflector_wall.h"
#include "corgi/corgi.h"
#include "corgi/tile.h"
#include "pybind11/numpy.h"
#include "thrust/device_vector.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace pic {

namespace mpi = mpi4cpp::mpi;

enum class ParticlePusher { boris, higuera_cary, faraday };
enum class FieldInterpolator { linear_1st, linear_1st_unrolled };
enum class CurrentDepositer { zigzag_1st, zigzag_1st_atomic };

struct ParticleStateBatch {
  using container_type = std::array<pybind11::array_t<double>, 3>;

  container_type pos;
  container_type vel;
};

/*! \brief pic tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from emf::Tile
 */

template<std::size_t D>
class Tile : virtual public emf::Tile<D> {
public:
  /// The type in which pos and vel are stored in.
  using value_type = ParticleContainer::value_type;

  static_assert(D == 3);

private:
  /// Tag which is shared between all the particles from this tile.
  std::size_t tile_tag_;
  // This should not be modified manually.
  std::map<std::size_t, std::size_t> next_particle_ordinals_;
  runko::prtc_id_type consume_next_id_(std::size_t ptype);

  std::map<std::size_t, ParticleContainer> particle_buffs_;
  ParticlePusher particle_pusher_;
  FieldInterpolator field_interpolator_;
  CurrentDepositer current_depositer_;

  std::vector<pic::reflector_wall> reflector_walls_ {};
  std::optional<runko::VecGrid<emf::YeeLattice::value_type>> reflector_correction_J_ {};
  std::map<std::size_t, std::size_t> amount_of_particles_to_be_send_;
  std::map<std::size_t, std::size_t> amount_of_particles_to_be_received_;


  /// particle type, direction -> span in subregion_particle_buff_
  ///
  /// These can not be std::span,
  /// because resizing subregion_particle_buff_ will invalidate pointers.
  std::map<
    std::size_t,
    std::map<runko::grid_neighbor<3>, std::pair<std::size_t, std::size_t>>>
    subregion_particle_spans_;
  thrust::device_vector<runko::ParticleState<value_type>> subregion_particle_buff_;

  void divide_particles_to_subregions();

  /// particle type -> incoming particles
  std::map<std::size_t, std::vector<std::span<const runko::ParticleState<value_type>>>>
    incoming_subregion_particles_ {};

public:
  /// Construct Tile based on the given config.
  ///
  /// FIXME: document how the initialization is done?
  ///
  /// In addition to emf::Tile ctor requirements,
  /// the given config has to contain values for:
  ///
  /// `qx`:    charge of x:th particle species (x is natural number)
  /// `mx`:    mass of x:th particle species (x is natural number)
  /// `c_omp`: simulation skin depth resolution
  /// `ppc`:   particles per cell per species
  ///
  /// `particle_pusher`:     scheme to update particles velocities and positions
  /// `fields_interpolator`: scheme to interpolate E and B fields to particles
  /// `current_depositer`:   scheme to depot current
  ///
  /// Note that particle charges and masses qx and mx are read in order: 0, 1, ...
  /// If a i:th mass and charge are missing, the search is stopped.
  explicit Tile(
    std::array<std::size_t, 3> tile_grid_idx,
    const toolbox::ConfigParser& config);

  // Has to be explicitly declared as a work around for hipcc bug.
  // see: https://github.com/llvm/llvm-project/issues/141592
  ~Tile() = default;

  std::array<std::vector<value_type>, 3> get_positions(std::size_t);
  std::array<std::vector<value_type>, 3> get_velocities(std::size_t);
  std::vector<runko::prtc_id_type> get_ids(std::size_t);


  using particle_generator =
    std::function<std::vector<runko::ParticleState<double>>(double, double, double)>;

  /// Inject particles based on given generator.
  ///
  /// Generator is called for each cell coordinates.
  ///
  /// Particle type is assumed to be configured.
  void inject_to_each_cell(std::size_t particle_type, particle_generator);

  /// Inject given particles.
  ///
  /// Particle type is assumed to be configured.
  void inject(std::size_t particle_type, std::vector<runko::ParticleState<double>>);

  using batch_array = pybind11::array_t<double>;
  using batch_particle_generator =
    std::function<ParticleStateBatch(batch_array, batch_array, batch_array)>;

  /// Inject particles based on given generator.
  ///
  /// Generator is called once with all cell coordinates.
  ///
  /// Particle type is assumed to be configured.
  void batch_inject_to_cells(std::size_t particle_type, batch_particle_generator);

  /// Inject particles in a stripe between x_left and x_right.
  ///
  /// Only cells whose x-coordinate falls within [x_left, x_right) are
  /// passed to the generator. Full y and z extent of the tile is used.
  /// No-op if the stripe does not overlap this tile.
  void batch_inject_in_x_stripe(
    std::size_t particle_type,
    batch_particle_generator pgen,
    double x_left,
    double x_right);

  /// Push particles updating their velocities and positions.
  void push_particles();

  /// Deposit current from all particls.
  void deposit_current();

  /// Sorts the particles in order to reduce cache misses.
  void sort_particles();

  /// Register a reflector wall on this tile.
  void register_reflector_wall(pic::reflector_wall wall);

  /// Reflect particles that crossed any registered reflector wall.
  ///
  /// Must be called after push_particles() and before deposit_current().
  /// Modifies particle positions/velocities for reflected particles
  /// and stores correction currents that deposit_current() will add.
  void reflect_particles();

  /// Update reflector wall locations by their velocity * cfl.
  void advance_reflector_walls();

  /// Interpolate E and B fields to given particle positions.
  ///
  /// Positions are in tile-local code units (same as ParticleContainer).
  /// Uses linear_1st interpolation.
  emf::YeeLattice::InterpolatedEB
    interpolate_fields_at(const runko::VecList<value_type>& positions) const;

  std::size_t number_of_species() const;

  std::size_t number_of_particles(std::size_t particle_type) const;

  /// Const access to particle container for a given species.
  const ParticleContainer& particles(std::size_t species) const;

  /// Returns average kinetic energy of given particle type.
  ///
  /// Kinetic energy is given in units of mc^2,
  /// where m is the mass corresponding to the particle type.
  ///
  /// Particle type is assumed to be configured.
  double total_kinetic_energy(std::size_t particle_type) const;

  std::vector<mpi4cpp::mpi::request>
    send_data(mpi4cpp::mpi::communicator& /*comm*/, int dest, int mode, int tag)
      override;

  std::vector<mpi4cpp::mpi::request>
    recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag)
      override;

  void local_communication_prelude(const int) override;
  void local_communication_postlude(const int) override;

  /// Get particles from haloregion of the other with comm_mode::pic_particle.
  ///
  /// Forward other communication modes to emf::Tile.
  /// Assumes that the other tile is pic::Tile or its descendant.
  void local_communication(
    const corgi::Tile<D>& /* other */,
    const std::array<int, D> dir_to_other,
    const int /* mode */
    ) override;
};


}  // namespace pic
