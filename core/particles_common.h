#pragma once

#include <array>

namespace runko {

enum class particle : std::size_t {
  electron = 0,
  ion      = 1,
  photon   = 2,
  proton   = 3,

};

using prtc_id_type                 = std::uint64_t;
static constexpr auto dead_prtc_id = std::numeric_limits<prtc_id_type>::max();

template <typename T>
struct ParticleState {
  using vec3 = std::array<T, 3>;

  vec3 pos;
  vec3 vel;
  prtc_id_type id;
};

static_assert(std::is_trivially_copyable_v<ParticleState<float>>);
}  // namespace runko
