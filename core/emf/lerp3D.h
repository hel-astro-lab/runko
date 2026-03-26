#pragma once

#include <type_traits>

namespace emf {

/// Linear interpolation of mds values to pos inside a unit cube.
///
/// Values in mds are assumed to be located at (i, j, k),
/// where i, j and k are either 0 or 1. Pos is assumed to be
/// inside the unit cube and measured from (0, 0, 0).
constexpr auto
  lerp3D(const auto& pos, const auto& mds)
{
  using MDS = std::remove_cvref_t<decltype(mds)>;
  static_assert(MDS::rank() == 3);
  static_assert(MDS::static_extent(0) == 2);
  static_assert(MDS::static_extent(1) == 2);
  static_assert(MDS::static_extent(2) == 2);

  constexpr auto lerp = [](const auto x, const auto A, const auto B) {
    return (1 - x) * A + x * B;
  };

  const auto c00 = lerp(pos[0], mds[0, 0, 0], mds[1, 0, 0]);
  const auto c01 = lerp(pos[0], mds[0, 0, 1], mds[1, 0, 1]);
  const auto c10 = lerp(pos[0], mds[0, 1, 0], mds[1, 1, 0]);
  const auto c11 = lerp(pos[0], mds[0, 1, 1], mds[1, 1, 1]);

  const auto c0 = lerp(pos[1], c00, c10);
  const auto c1 = lerp(pos[1], c01, c11);

  return lerp(pos[2], c0, c1);
}

}  // namespace emf
