#pragma once

#include "tyvi/mdspan.h"

#include <algorithm>
#include <compare>
#include <concepts>
#include <cstddef>
#include <optional>
#include <print>
#include <stdexcept>
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
  emf_J_exchange
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
  std::is_enum_v<E> and std::convertible_to<T, std::underlying_type_t<E>>;

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

template<std::size_t rank>
struct [[nodiscard]] grid_neighbor {
  using value_type = std::int8_t;
  std::array<value_type, rank> direction;

  constexpr std::strong_ordering operator<=>(const grid_neighbor&) const = default;

  template<std::convertible_to<value_type>... T>
  explicit constexpr grid_neighbor(T&&... args) :
    direction { static_cast<value_type>(std::forward<T>(args))... }
  {
    const auto ok = std::ranges::all_of(direction, [](const auto x) {
      const auto a = x == value_type { -1 };
      const auto b = x == value_type { 0 };
      const auto c = x == value_type { 1 };

      return a or b or c;
    });

    if(not ok) {
      throw std::logic_error {
        "All components in grid_neighbor has to be in {-1, 0, 1}."
      };
    }
  }

  template<std::convertible_to<value_type> T>
  constexpr grid_neighbor(const std::array<T, rank>& arr)
  {
    *this = [&]<std::size_t... I>(std::index_sequence<I...>) {
      return grid_neighbor(arr[I]...);
    }(std::make_index_sequence<rank>());
  }


  constexpr grid_neighbor inverted() const
  {
    return [this]<std::size_t... I>(std::index_sequence<I...>) {
      return grid_neighbor(-direction[I]...);
    }(std::make_index_sequence<rank>());
  }

  /* nieghbor_index one-to-one maps directions to numbers: 0, 1, ..., 3^rank - 1 */

  static constexpr auto mapping = []<std::size_t... I>(std::index_sequence<I...>) {
    // This probably could have value_type as index type but due to
    // https://github.com/hel-astro-lab/tyvi/issues/11
    // we have to use wider type.
    using E = std::extents<unsigned int, (I * 0 + 3)...>;
    return std::layout_right::mapping { E {} };
  }(std::make_index_sequence<rank>());

  static constexpr auto inverse_mapping = tyvi::sstd::index_space_view(mapping);

  constexpr std::uint8_t neighbor_index() const
  {
    // ugly syntax until template for in C++26 :(
    const auto index = [&]<std::size_t... I>(std::index_sequence<I...>) {
      return mapping((direction[I] + 1)...);
    }(std::make_index_sequence<rank>());
    static_assert(
      rank <= 5uz,
      "rank-6 or larger neighbor_index does not fit into std::uint8_t");
    return static_cast<std::uint8_t>(index);
  }

  static constexpr grid_neighbor from_index(const std::uint8_t index)
  {
    // ugly syntax until template for in C++26 :(
    return [&]<std::size_t... I>(std::index_sequence<I...>) {
      const auto dir = inverse_mapping[index];
      return grid_neighbor((dir[I] - 1)...);
    }(std::make_index_sequence<rank>());
  }
};

static_assert(
  grid_neighbor<2>(-1, -1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(-1, -1).neighbor_index()));
static_assert(
  grid_neighbor<2>(-1, 0) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(-1, 0).neighbor_index()));
static_assert(
  grid_neighbor<2>(-1, 1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(-1, 1).neighbor_index()));
static_assert(
  grid_neighbor<2>(0, -1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(0, -1).neighbor_index()));
static_assert(
  grid_neighbor<2>(0, 0) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(0, 0).neighbor_index()));
static_assert(
  grid_neighbor<2>(0, 1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(0, 1).neighbor_index()));
static_assert(
  grid_neighbor<2>(1, -1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(1, -1).neighbor_index()));
static_assert(
  grid_neighbor<2>(1, 0) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(1, 0).neighbor_index()));
static_assert(
  grid_neighbor<2>(1, 1) ==
  grid_neighbor<2>::from_index(grid_neighbor<2>(1, 1).neighbor_index()));

}  // namespace runko
