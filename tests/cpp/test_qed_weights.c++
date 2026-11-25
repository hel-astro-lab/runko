// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/weights.h"
#include "runko/tools/random.h"

#include <boost/ut.hpp>
#include <cmath>

namespace {

using namespace boost::ut;

constexpr bool
  close(const double a, const double b, const double rtol)
{ return std::abs(a - b) <= rtol * std::abs(b); }

[[maybe_unused]]
const suite<"qed weight bookkeeping"> _ = [] {
  "event size times pair weight is the product of the weights"_test = [] {
    auto rng = runko::Rng { runko::rng_key(1, 0, 0, 0), 0 };

    for(auto n = 0; n < 100000; ++n) {
      const auto wi = std::exp(10.0 * rng.uniform() - 5.0);
      const auto wj = std::exp(10.0 * rng.uniform() - 5.0);
      const auto fk = n % 2 == 0 ? 1.0 : std::exp(8.0 * rng.uniform());
      const auto S  = qed::event_size(wi, wj, fk);

      // probability times size is the reaction count of the pair
      expect(close(S * qed::pair_weight(wi, wj, fk), wi * wj, 1e-12));

      // the incident never gives up more than it has; without merging the smaller
      // weight moves
      expect(S <= wi);
      if(fk == 1.0) { expect(S == std::min(wi, wj)); }

      // the partner takes at most fk exchanges
      expect(S / wj <= fk * (1.0 + 1e-12));
    }
  };

  "fkin is one without kappa and never below"_test = [] {
    expect(qed::fkin(0.0, 10.0, 1e-9) == 1.0);
    expect(qed::fkin(0.1, 1.0, 1e-3) == 1.0);
    expect(qed::fkin(0.1, 3.0, 10.0) == 1.0);
    expect(close(qed::fkin(0.1, 11.0, 1e-3), 1e3, 1e-12));
  };

  "copies conserve the weight, each draw above one and on average below"_test = [] {
    auto rng = runko::Rng { runko::rng_key(2, 0, 0, 0), 0 };

    for(const auto n: { 0.03, 0.4, 1.0, 2.6, 17.5, 80.0 }) {
      constexpr auto s = 0.7, draws = 200000.0;
      constexpr auto n_cap = 32;
      auto sum_n = 0.0, sum_w = 0.0;

      for(auto k = 0; k < static_cast<int>(draws); ++k) {
        const auto cp = qed::copies(n, s, rng, n_cap);
        if(n >= 1.0) { expect(close(cp.n_r * cp.w_each, s, 1e-12)); }
        sum_n += cp.n_r;
        sum_w += cp.n_r * cp.w_each;
      }

      const auto nc = std::min(n, static_cast<double>(n_cap));
      expect(close(sum_n / draws, nc, 4.0 / std::sqrt(nc * draws) + 1e-12));
      expect(close(sum_w / draws, s, 4.0 / std::sqrt(nc * draws) + 1e-12));
    }
  };

  "the merged emission stays a small recoil"_test = [] {
    for(const auto gam: { 3.0, 30.0, 3e3 }) {
      for(const auto xm: { 1e-6, 1e-2, 1.0, 50.0 }) {
        const auto f = qed::facc_syn(gam, xm, 3.0, 0.1);
        expect(f >= 1.0);
        expect(f == 1.0 or f * xm <= 0.1 * (gam - 1.0) * (1.0 + 1e-12));
        expect(qed::facc_syn(gam, xm, 0.0, 0.1) == 1.0);
      }
    }
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{ return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv })); }
