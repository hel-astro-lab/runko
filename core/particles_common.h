#pragma once

#include "tools/vector.h"

#include <array>
#include <concepts>
#include <type_traits>

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

template<typename F, typename T>
concept EB_interpolator =
  (std::is_trivially_copyable_v<F>) and std::regular_invocable<F, toolbox::Vec3<T>> and
  requires(std::invoke_result_t<F, toolbox::Vec3<T>> eb) {
    { eb.E } -> std::convertible_to<toolbox::Vec3<T>>;
    { eb.B } -> std::convertible_to<toolbox::Vec3<T>>;
  };

}  // namespace runko
