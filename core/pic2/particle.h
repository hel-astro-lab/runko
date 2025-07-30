#pragma once

#include "core/emf2/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "core/particles_common.h"
#include "thrust/device_vector.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <ranges>

namespace pic2 {

struct ParticleContainerArgs {
  std::size_t N;
  double charge, mass;
};

class [[nodiscard]] ParticleContainer {
public:
  /// The type in which pos and vel are stored in.
  using value_type = float;

private:
  using E = std::dextents<std::size_t, 1>;

  runko::VecList<value_type> pos_;
  runko::VecList<value_type> vel_;

  double charge_;
  double mass_;

public:
  explicit ParticleContainer(ParticleContainerArgs);

  /// Returns the number of particles.
  std::size_t size() const;
  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();

  [[nodiscard]]
  auto span_pos(this auto& self)
  {
    return self.pos_.span();
  }

  [[nodiscard]]
  auto span_vel(this auto& self)
  {
    return self.vel_.span();
  }

  /// Splits container to 27 subcontainers along given divider lines.
  ///
  /// Container pointed by this is replaced by the middle container
  /// and the rest are returned with corresponding direction.
  ///
  /// global_{mins,maxs} are required for global periodic boundary condition.
  std::vector<std::pair<std::array<int, 3>, ParticleContainer>> split_to_subregions(
    std::array<value_type, 2> x_dividers,
    std::array<value_type, 2> y_dividers,
    std::array<value_type, 2> z_dividers,
    std::array<value_type, 3> global_mins,
    std::array<value_type, 3> global_maxs);

  /// Add particles from other container.
  void add_particles(const ParticleContainer& other);

  /// Way to initialize particles. Known to be slow.
  template<std::ranges::forward_range R>
    requires std::
      convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState>
    void add_particles(const R& new_particles)
  {
    const auto N = static_cast<std::size_t>(std::ranges::distance(new_particles));

    auto args  = this->args();
    args.N     = N;
    auto added = ParticleContainer(args);

    const auto added_pos_smds = added.pos_.staging_mds();
    const auto added_vel_smds = added.vel_.staging_mds();

    for(const auto [i, p]: std::views::enumerate(new_particles)) {
      for(const auto j: std::views::iota(0uz, 3uz)) {
        added_pos_smds[i][j] = p.pos[j];
        added_vel_smds[i][j] = p.vel[j];
      }
    }

    auto w1 = tyvi::mdgrid_work {}.sync_from_staging(added.pos_);
    auto w2 = tyvi::mdgrid_work {}.sync_from_staging(added.vel_);
    tyvi::when_all(w1, w2).wait();

    this->add_particles(added);
  }

  using InterpolatedEB_function =
    std::function<emf2::YeeLattice::InterpolatedEB(const runko::VecList<value_type>&)>;

  /// Push particles velocities and positions using boris scheme.
  void push_particles_boris(double cfl, const InterpolatedEB_function&);

  emf2::YeeLattice::CurrentContributions current_zigzag_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const;
};


}  // namespace pic2
