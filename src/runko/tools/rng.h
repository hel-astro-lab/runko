// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Counter-based random numbers for use inside tyvi for_each kernels.
// Stateless and branch-free: a (counter, key) pair maps to the same output on every
// backend, lane and MPI rank, so nothing is seeded or stored per particle.
// Integer ops only: 32x32->64 multiplies, xors and adds.

#include "runko/tools/math.h"

#include <array>
#include <cstdint>
#include <numbers>

namespace toolbox {

/// Philox4x32-10 (Salmon+11): 4 output words from a 128-bit counter and 64-bit key.
constexpr std::array<std::uint32_t, 4>
  philox4x32(std::array<std::uint32_t, 4> ctr, std::array<std::uint32_t, 2> key)
{
  constexpr std::uint32_t M0 = 0xD2511F53u, M1 = 0xCD9E8D57u;  // round multipliers
  constexpr std::uint32_t W0 = 0x9E3779B9u, W1 = 0xBB67AE85u;  // Weyl key bumps
  for(int r = 0; r < 10; ++r) {
    const std::uint64_t p0 = std::uint64_t { M0 } * ctr[0];
    const std::uint64_t p1 = std::uint64_t { M1 } * ctr[2];
    ctr = { static_cast<std::uint32_t>(p1 >> 32) ^ ctr[1] ^ key[0],
            static_cast<std::uint32_t>(p1),
            static_cast<std::uint32_t>(p0 >> 32) ^ ctr[3] ^ key[1],
            static_cast<std::uint32_t>(p0) };
    key[0] += W0;
    key[1] += W1;
  }
  return ctr;
}

// Random123 known-answer vector
static_assert(philox4x32({ 0, 0, 0, 0 }, { 0, 0 })
              == std::array<std::uint32_t, 4> { 0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u });

/// Uniform in (0, 1]; never 0 so log() is safe.
template<typename vt>
constexpr vt
  uniform01(const std::uint32_t r)
{
  return (static_cast<vt>(r) + vt { 1 }) * vt { 2.3283064365386963e-10 };  // 2^-32
}

/// Four standard normals for particle `id` at counter `cntr` of species `sp` in run `seed`.
template<typename vt>
constexpr std::array<vt, 4>
  normal4(
    const std::uint64_t id,
    const std::uint64_t cntr,
    const std::uint32_t sp,
    const std::uint32_t seed)
{
  const auto r = philox4x32(
    { static_cast<std::uint32_t>(id),
      static_cast<std::uint32_t>(id >> 32),
      static_cast<std::uint32_t>(cntr),
      static_cast<std::uint32_t>(cntr >> 32) },
    { sp, seed });
  // Box-Muller: 2 log + 2 sincos; Irwin-Hall sum of 12 uniforms if this shows in profiles
  constexpr vt twopi = vt { 2 } * std::numbers::pi_v<vt>;
  const vt m0 = sstd::sqrt(vt { -2 } * sstd::log(uniform01<vt>(r[0])));
  const vt a0 = twopi * uniform01<vt>(r[1]);
  const vt m1 = sstd::sqrt(vt { -2 } * sstd::log(uniform01<vt>(r[2])));
  const vt a1 = twopi * uniform01<vt>(r[3]);
  return { m0 * sstd::cos(a0), m0 * sstd::sin(a0), m1 * sstd::cos(a1), m1 * sstd::sin(a1) };
}

}  // namespace toolbox
