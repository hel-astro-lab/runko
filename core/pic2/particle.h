#pragma once

#include "tyvi/mdgrid.h"

#include <cstddef>

namespace pic2 {

struct ParticleContainerArgs {
  std::size_t N;
  double charge, mass;
};

template<std::size_t D>
class [[nodiscard]] ParticleContainer {

  using E = std::dextents<std::size_t, 1>;
  static constexpr auto ScalarElement =
    tyvi::mdgrid_element_descriptor<float> { .rank = 0, .dim = 3 };
  static constexpr auto VecElement =
    tyvi::mdgrid_element_descriptor<float> { .rank = 1, .dim = 3 };

  using VecGrid    = tyvi::mdgrid<VecElement, E>;
  using ScalarGrid = tyvi::mdgrid<VecElement, E>;

  VecGrid pos_;
  VecGrid vel_;
  ScalarGrid weights_;

  double charge_;
  double mass_;

public:
  explicit ParticleContainer(const ParticleContainerArgs args) :
    pos_(args.N),
    vel_(args.N),
    weights_(args.N),
    charge_ { args.charge },
    mass_ { args.mass }
  {
  }
};


}  // namespace pic2
