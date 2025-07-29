#pragma once

#include <concepts>
#include <optional>
#include <type_traits>
#include <utility>

/* Corgi {send,recv}_data works with different modes.
   These modes are specified with `int`,
   which is then passed to {send,recv}_data members of tiles.

   This header contains enum declerations for each mode.
   Additionally there are python bindings for these in pytools module,
   so that python and C++ can have common names for different modes. */

namespace runko {
enum class comm_mode : int {
  emf_E              = 1,
  emf_B              = 2,
  emf_J              = 0,
  pic_particle       = 3,
  pic_particle_extra = 4,
  number_of_particles,
};


/// Returns a communication which has to be done before given communication.
[[nodiscard]] constexpr std::optional<comm_mode>
  virtual_tile_sync_handshake_mode(const comm_mode mode)
{
  switch(mode) {
    case comm_mode::pic_particle: return comm_mode::number_of_particles;
    default: return {};
  }
}


template<typename T, typename E>
concept equality_comparable_with_enum =
  (std::is_enum_v<E>) and std::equality_comparable_with<T, std::underlying_type_t<E>>;

template<typename Rhs>
  requires(not std::same_as<Rhs, comm_mode>) and
          equality_comparable_with_enum<Rhs, comm_mode>
[[nodiscard]]
constexpr bool
  operator==(const comm_mode lhs, const Rhs rhs)
{
  return std::to_underlying(lhs) == rhs;
}

template<typename Lhs>
  requires(not std::same_as<Lhs, comm_mode>) and
          equality_comparable_with_enum<Lhs, comm_mode>
[[nodiscard]]
constexpr bool
  operator==(const Lhs lhs, const comm_mode rhs)
{
  return lhs == std::to_underlying(rhs);
}

static_assert(comm_mode::emf_E == std::to_underlying(comm_mode::emf_E));
static_assert(std::to_underlying(comm_mode::emf_E) == comm_mode::emf_E);


}  // namespace runko
