// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/ut.hpp>

namespace {

constexpr auto
  sum(auto... values)
{
  return (values + ...);
}

auto
  call_expect_with(const bool x)
{
  boost::ut::expect(x);
}

using namespace boost::ut;
[[maybe_unused]]
const suite<"unit testing"> _ = [] {
  "sum"_test = [] {
    expect(sum(0) == 0_i);
    expect(sum(1, 2) == 3_i);
    expect(sum(1, 2) > 0_i and 4_i == sum(2, 2));
  };

  "subfunctions"_test = [] { call_expect_with(true); };
};

}  // namespace

int
  main(int argc, const char** argv)
{
  return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv }));
}
