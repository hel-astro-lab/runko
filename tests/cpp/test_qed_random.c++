// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "thrust/device_vector.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

using namespace boost::ut;

[[maybe_unused]]
const suite<"counter-based rng"> _ = [] {
  "uniform moments and bins"_test = [] {
    const auto key   = runko::rng_key(1, 2, 0, 3);
    auto rng         = runko::Rng { key, 0 };
    constexpr auto n = 1000000;
    auto sum = 0.0, sum2 = 0.0, corr = 0.0, prev = 0.5;
    auto bins = std::array<int, 16> {};

    for(auto i = 0; i < n; ++i) {
      const auto u = rng.uniform();
      expect(u > 0.0 and u < 1.0);
      sum += u;
      sum2 += u * u;
      corr += (u - 0.5) * (prev - 0.5);
      prev = u;
      ++bins[static_cast<std::size_t>(u * 16.0)];
    }

    const auto mean = sum / n;
    const auto var  = sum2 / n - mean * mean;
    expect(std::abs(mean - 0.5) < 2e-3);
    expect(std::abs(var - 1.0 / 12.0) < 2e-3);
    expect(std::abs(corr / n) < 2e-3);
    for(const auto b: bins) { expect(std::abs(b - n / 16.0) < 0.03 * n / 16.0); }
  };

  "a draw depends on (key, lane, index) only"_test = [] {
    const auto key = runko::rng_key(7, 11, 1, 2);
    auto a         = runko::Rng { key, 5 };
    auto b         = runko::Rng { key, 5 };
    for(auto i = 0; i < 100; ++i) { expect(a.uniform() == b.uniform()); }

    auto c = runko::Rng { key, 6 };
    auto d = runko::Rng { runko::rng_key(7, 11, 1, 3), 5 };
    auto e = runko::Rng { key, 5 };
    expect(c.uniform() != e.uniform());
    expect(d.uniform() != runko::Rng { key, 5 }.uniform());
  };

  "lanes are independent"_test = [] {
    // the first draws of consecutive lanes look like one uniform sample
    const auto key   = runko::rng_key(3, 0, 0, 0);
    constexpr auto n = 200000;
    auto sum = 0.0, corr = 0.0, prev = 0.5;

    for(auto lane = 0; lane < n; ++lane) {
      const auto u = runko::Rng { key, static_cast<std::uint64_t>(lane) }.uniform();
      sum += u;
      corr += (u - 0.5) * (prev - 0.5);
      prev = u;
    }

    expect(std::abs(sum / n - 0.5) < 3e-3);
    expect(std::abs(corr / n) < 3e-3);
  };

  "normal moments"_test = [] {
    auto rng         = runko::Rng { runko::rng_key(5, 0, 0, 0), 0 };
    constexpr auto n = 400000;
    auto s = 0.0, s2 = 0.0, s4 = 0.0;

    for(auto i = 0; i < n; ++i) {
      const auto z = rng.normal();
      s += z;
      s2 += z * z;
      s4 += z * z * z * z;
    }

    expect(std::abs(s / n) < 5e-3);
    expect(std::abs(s2 / n - 1.0) < 1e-2);
    expect(std::abs(s4 / n - 3.0) < 6e-2);
  };

  "integer range"_test = [] {
    auto rng = runko::Rng { runko::rng_key(9, 0, 0, 0), 0 };
    auto cnt = std::array<int, 7> {};
    for(auto i = 0; i < 70000; ++i) { ++cnt[rng.integer(7)]; }
    for(const auto c: cnt) { expect(std::abs(c - 10000) < 400); }
  };
};

[[maybe_unused]]
const suite<"portable atomics"> atomics = [] {
  "compare_exchange: of many lanes racing on one slot exactly one wins"_test = [] {
    // the winner stored its value and saw true; every loser saw false and left it
    constexpr auto n = 100000;
    auto slot        = thrust::device_vector<std::uint64_t>(1, 0);
    auto won         = thrust::device_vector<int>(n, 0);
    const auto s     = thrust::raw_pointer_cast(slot.data());
    const auto wp    = thrust::raw_pointer_cast(won.data());
    const auto w     = tyvi::mdgrid_work {};
    const auto i0    = thrust::counting_iterator<int>(0);

    thrust::for_each(w.on_this(), i0, i0 + n, [=](const int i) {
      wp[i] =
        sstd::atomic_compare_exchange(s, std::uint64_t { 0 }, std::uint64_t(i) + 1) ? 1
                                                                                    : 0;
    });
    w.wait();

    const auto h      = thrust::host_vector<int>(won);
    const auto winner = std::find(h.begin(), h.end(), 1) - h.begin();
    expect(std::count(h.begin(), h.end(), 1) == 1);
    expect(std::uint64_t { slot[0] } == static_cast<std::uint64_t>(winner) + 1);
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{
  return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv }));
}
