#pragma once

#include "core/mdgrid_common.h"
#include "core/particles_common.h"

#include <concepts>
#include <cstddef>
#include <ranges>

namespace pic2 {

struct ParticleContainerArgs {
  std::size_t N;
  double charge, mass;
};

class [[nodiscard]] ParticleContainer {
public:
  /// The type in which pos, vel and weights are stored in.
  using value_type = float;

private:
  using E = std::dextents<std::size_t, 1>;

  runko::VecList<float> pos_;
  runko::VecList<float> vel_;
  runko::ScalarList<float> weights_;

  double charge_;
  double mass_;

public:
  explicit ParticleContainer(ParticleContainerArgs);

  /// Returns the number of particles.
  std::size_t size() const;
  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();
  std::vector<value_type> get_weights();

  /// Add particles from other container.
  ///
  /// Returns tyvi:::mdgrid_work representing the ongoing work.
  void add_particles(ParticleContainer& other);


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

    const auto added_pos_smds     = added.pos_.staging_mds();
    const auto added_vel_smds     = added.vel_.staging_mds();
    const auto added_weights_smds = added.weights_.staging_mds();

    for(const auto [i, p]: std::views::enumerate(new_particles)) {
      for(const auto j: std::views::iota(0uz, 3uz)) {
        added_pos_smds[i][j] = p.pos[j];
        added_vel_smds[i][j] = p.vel[j];
      }
      added_weights_smds[i][] = p.weight;
    }

    auto w1 = tyvi::mdgrid_work {}.sync_from_staging(added.pos_);
    auto w2 = tyvi::mdgrid_work {}.sync_from_staging(added.vel_);
    auto w3 = tyvi::mdgrid_work {}.sync_from_staging(added.weights_);
    tyvi::when_all(w1, w2, w3).wait();

    this->add_particles(added);
  }
};


}  // namespace pic2
