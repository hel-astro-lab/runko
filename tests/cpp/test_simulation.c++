// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/simulation.h"
#include "tyvi/mdspan.h"
#include "tyvi/sstd.h"

#include <boost/ut.hpp>
#include <ranges>
#include <set>

namespace {
using namespace boost::ut;

const auto s = [] {
  "initializing cartesian tiles"_test = [] {
    auto sim = runko::simulation_context {};

    runko::add_cartesian_tiles<3>(sim, { 3uz, 2uz, 4uz });

    using arr = decltype(runko::cartesian_index<3>::data);

    auto seen = std::set<arr> {};

    for(const auto [_, idx]: sim.view_tiles<const runko::cartesian_index<3>>()) {
      seen.insert(idx.data);
    }

    for(const auto idx: tyvi::sstd::index_space_view(
          std::layout_right::mapping(std::dextents<std::size_t, 3>(3, 2, 4)))) {
      expect(seen.contains(idx));
    }
  };
};
}  // namespace


int
  main(int argc, const char** argv)
{
  [[maybe_unused]]
  const suite<"simulation"> _ = s;
  return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv }));
}
