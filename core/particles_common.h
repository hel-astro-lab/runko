#pragma once

#include <array>

namespace runko {

enum class particle : std::size_t {
  electron = 0,
  ion      = 1,
  photon   = 2,

};

struct ParticleState {
  using vec3 = std::array<double, 3>;

  vec3 pos;
  vec3 vel;
  double weight;
};

}  // namespace runko
