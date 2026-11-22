// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/tools/math.h"

#include <cstdint>

namespace runko {

// Counter-based random numbers: a draw is a hash of (stream key, counter) and nothing
// else, so a kernel lane owns a stream fixed by (seed, step, sub-step, kernel, lane)
// and no state is carried between lanes or launches. Reproducible whatever the launch
// order, and free of the per-lane state a Mersenne twister would need.

/// splitmix64 finalizer: a bijection of the 64-bit integers with full avalanche.
[[nodiscard]]
constexpr std::uint64_t
  mix64(std::uint64_t z)
{
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

/// Stream key of one kernel launch: seed, time step, sub-step and kernel id folded
/// together.
[[nodiscard]]
constexpr std::uint64_t
  rng_key(
    const std::uint64_t seed,
    const std::uint64_t step,
    const std::uint64_t sub,
    const std::uint64_t kernel)
{ return mix64(mix64(mix64(seed) ^ step) ^ (sub << 32) ^ kernel); }

/// The stream of one lane: uniform, normal and integer draws numbered by a counter.
struct [[nodiscard]] Rng {
  std::uint64_t key;  // the stream of the launch (seed, step, sub-step, kernel)
  std::uint64_t ctr;  // the next draw of this lane

  // counter space of one lane
  static constexpr std::uint64_t draws_per_lane = std::uint64_t { 1 } << 24;

  constexpr Rng(const std::uint64_t stream_key, const std::uint64_t lane) :
    key { stream_key },
    ctr { lane * draws_per_lane }
  {
  }

  /// Next raw 64-bit word of the stream.
  [[nodiscard]]
  constexpr std::uint64_t next()
  {
    const auto x = mix64(ctr++ + key);  // first mix of counter and key
    return mix64(
      x ^ ((key << 32) | (key >> 32)));  // second mix against the rotated key
  }

  /// Uniform on the open interval (0, 1): 53 random bits plus a half-ulp offset, never
  /// 0 or 1.
  [[nodiscard]]
  constexpr double uniform()
  { return (static_cast<double>(next() >> 11) + 0.5) * 0x1.0p-53; }

  /// Uniform integer in [0, n).
  [[nodiscard]]
  constexpr std::uint64_t integer(const std::uint64_t n)
  { return static_cast<std::uint64_t>(uniform() * static_cast<double>(n)) % n; }

  /// Standard normal deviate by the Box-Muller transform.
  [[nodiscard]]
  constexpr double normal()
  {
    const auto r   = sstd::sqrt(-2.0 * sstd::log(uniform()));  // Rayleigh radius
    const auto phi = sstd::TWOPI * uniform();                  // uniform angle
    return r * sstd::cos(phi);
  }
};

}  // namespace runko
