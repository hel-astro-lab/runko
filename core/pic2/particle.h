#pragma once

#include "tyvi/mdgrid.h"

#include <cstddef>

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
  static constexpr auto scalar_element =
    tyvi::mdgrid_element_descriptor<value_type> { .rank = 0, .dim = 3 };
  static constexpr auto vec_element =
    tyvi::mdgrid_element_descriptor<value_type> { .rank = 1, .dim = 3 };

  using ScalarGrid = tyvi::mdgrid<scalar_element, E>;
  using VecGrid    = tyvi::mdgrid<vec_element, E>;

  VecGrid pos_;
  VecGrid vel_;
  ScalarGrid weights_;

  double charge_;
  double mass_;

public:
  ParticleContainer 
  explicit ParticleContainer(ParticleContainerArgs);

  /// Returns the number of particles.
  std::size_t size() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();
  std::vector<value_type> get_weights();
};


}  // namespace pic2
