// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/ut.hpp>

#include <cmath>

#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed_kinematics.h"
#include "runko/tools/random.h"

namespace {

using namespace boost::ut;
using qed::Vec3;

constexpr bool
  close(const double a, const double b, const double rtol, const double atol = 0.0)
{
  return std::abs(a - b) <= rtol * std::abs(b) + atol;
}

Vec3
  random_u(runko::Rng& rng, const double pmax)
{
  const auto p = pmax * std::pow(rng.uniform(), 0.5);
  return p * qed::iso_dir(rng);
}

double
  gam(const Vec3& u)
{
  return std::sqrt(1.0 + toolbox::dot(u, u));
}

[[maybe_unused]]
const suite<"two-body qed physics"> _ = [] {
  "pinned values from the python modules"_test = [] {
    expect(close(qed::comp::sigma_kn(1e-4), 0.9998000519867, 1e-11));
    expect(close(qed::comp::sigma_kn(5e-3), 0.9901283375, 1e-12));
    expect(close(qed::comp::sigma_kn(0.02), 0.9619785920554768, 1e-9));  // closed form cancels here
    expect(close(qed::comp::sigma_kn(0.5), 0.5627857603804511, 1e-12));
    expect(close(qed::comp::sigma_kn(2.0), 0.3141151957093031, 1e-12));
    expect(close(qed::comp::sigma_kn(50.0), 0.03755080938713749, 1e-12));
    expect(close(qed::ann::bracket(1.0), 2.0000000000000187, 1e-12));
    expect(close(qed::ann::bracket(1.001), 2.0019988007233906, 1e-12));
    expect(close(qed::ann::bracket(1.5), 2.7700523262672596, 1e-12));
    expect(close(qed::ann::bracket(3.0), 3.855477641542537, 1e-12));
    expect(close(qed::ann::bracket(20.0), 6.267961140237268, 1e-12));
    expect(close(qed::ann::bracket(1e4), 17.812516935563163, 1e-12));
    expect(qed::gg::sigma_gg(1.5) == 0.0);
    expect(qed::gg::sigma_gg(2.0) == 0.0);
    expect(close(qed::gg::sigma_gg(2.5), 0.18582075908295057, 1e-12));
    expect(close(qed::gg::sigma_gg(3.94), 0.2556393223638441, 1e-12));
    expect(close(qed::gg::sigma_gg(10.0), 0.17502503191665705, 1e-12));
    expect(close(qed::gg::sigma_gg(1e3), 0.004960566610367219, 1e-12));
    expect(close(
      qed::comp::rate(Vec3(0.7, 0.1, -0.4), Vec3(0.3, -0.2, 0.5)),
      0.22095192969772925,
      1e-12));
    expect(close(
      qed::ann::rate(Vec3(1.2, 0.3, -0.4), Vec3(-0.5, 0.9, 0.2)),
      0.1480432202450246,
      1e-12));
    expect(close(
      qed::gg::rate(Vec3(1.5, 0.2, -0.3), Vec3(-1.1, 0.8, 0.6)),
      0.22088514163908554,
      1e-12));
    expect(qed::gg::rate(Vec3(1.5, 0.2, -0.3), Vec3(0.2, 0.1, 0.05)) == 0.0);
  };

  "suprema respect SHMAX"_test = [] {
    auto rng = runko::Rng { runko::rng_key(1, 0, 0, 0), 0 };
    for(auto n = 0; n < 200000; ++n) {
      const auto k = random_u(rng, 20.0), u = random_u(rng, 20.0);
      expect(qed::comp::rate(k, u) <= qed::comp::SHMAX);
      expect(qed::ann::rate(k, u) <= qed::ann::SHMAX);
      expect(qed::gg::rate(k, u) <= qed::gg::SHMAX);
    }
  };

  "boost inverts and preserves the invariant mass"_test = [] {
    auto rng = runko::Rng { runko::rng_key(2, 0, 0, 0), 0 };
    for(auto n = 0; n < 1000; ++n) {
      const auto u  = random_u(rng, 5.0);
      const auto p  = random_u(rng, 3.0);
      const auto e  = std::sqrt(0.25 + toolbox::dot(p, p));
      const auto b  = u / gam(u);
      const auto A  = qed::boost(b, e, p);
      const auto B  = qed::boost((-1.0) * b, A.e, A.p);
      expect(close(B.e, e, 1e-12));
      for(auto i = 0uz; i < 3uz; ++i) { expect(close(B.p[i], p[i], 1e-12, 1e-14)); }
      expect(close(A.e * A.e - toolbox::dot(A.p, A.p), 0.25, 1e-10, 1e-12));
    }
  };

  "perp_basis is orthonormal, also for a parallel reference"_test = [] {
    auto rng = runko::Rng { runko::rng_key(3, 0, 0, 0), 0 };
    for(auto n = 0; n < 100; ++n) {
      const auto kv = qed::iso_dir(rng);
      for(const auto& ref: { qed::iso_dir(rng), kv, (-1.0) * kv }) {
        const auto t = qed::perp_basis(kv, ref);
        expect(close(toolbox::norm(t.i), 1.0, 1e-12));
        expect(close(toolbox::norm(t.j), 1.0, 1e-12));
        expect(close(toolbox::dot(t.i, t.j), 0.0, 0.0, 1e-12));
        expect(close(toolbox::dot(t.i, t.k), 0.0, 0.0, 1e-12));
        expect(close(toolbox::dot(t.j, t.k), 0.0, 0.0, 1e-12));
        const auto c = toolbox::cross(t.i, t.j);
        for(auto i = 0uz; i < 3uz; ++i) { expect(close(c[i], t.k[i], 0.0, 1e-12)); }
      }
    }
  };

  "compton sampler conserves four-momentum, few trips"_test = [] {
    auto rng   = runko::Rng { runko::rng_key(4, 0, 0, 0), 0 };
    auto trips = 0.0;
    auto worst = 0;
    constexpr auto n = 20000;
    for(auto m = 0; m < n; ++m) {
      const auto k = std::exp(-6.0 + 8.0 * rng.uniform()) * qed::iso_dir(rng);
      const auto u = random_u(rng, 30.0);
      auto ntry    = 0;
      const auto P = qed::comp::sample(k, u, rng, ntry);
      trips += ntry;
      worst = std::max(worst, ntry);
      const auto e0 = toolbox::norm(k) + gam(u);
      const auto e1 = toolbox::norm(P.k1) + gam(P.u1);
      expect(close(e1, e0, 1e-9)) << "energy";
      const auto p0 = k + u, p1 = P.k1 + P.u1;
      for(auto i = 0uz; i < 3uz; ++i) { expect(close(p1[i], p0[i], 1e-9, 1e-12)); }
    }
    expect(trips / n < 2.5);
    expect(worst < 40);
  };

  "compton reduces to Thomson: <mu^2> = 2/5 for a soft photon on a lepton at rest"_test =
    [] {
      auto rng = runko::Rng { runko::rng_key(5, 0, 0, 0), 0 };
      auto s2  = 0.0;
      constexpr auto n = 200000;
      const auto k = Vec3(1e-4, 0.0, 0.0);
      for(auto m = 0; m < n; ++m) {
        auto ntry    = 0;
        const auto P = qed::comp::sample(k, Vec3(0.0, 0.0, 0.0), rng, ntry);
        const auto mu = P.k1[0] / toolbox::norm(P.k1);
        s2 += mu * mu;
      }
      expect(close(s2 / n, 0.4, 0.01));
    };

  "annihilation sampler conserves four-momentum"_test = [] {
    auto rng   = runko::Rng { runko::rng_key(6, 0, 0, 0), 0 };
    auto trips = 0.0;
    constexpr auto n = 20000;
    for(auto m = 0; m < n; ++m) {
      const auto um = random_u(rng, 20.0), up = random_u(rng, 20.0);
      auto ntry     = 0;
      const auto P  = qed::ann::sample(um, up, rng, ntry);
      trips += ntry;
      expect(close(toolbox::norm(P.k0) + toolbox::norm(P.k1), gam(um) + gam(up), 1e-9));
      const auto p0 = um + up, p1 = P.k0 + P.k1;
      for(auto i = 0uz; i < 3uz; ++i) { expect(close(p1[i], p0[i], 1e-9, 1e-12)); }
    }
    expect(trips / n < 1.3);
    // a pair with no net momentum: two photons of energy gamma each, back to back, and
    // for the pair exactly at rest an isotropic direction
    auto c2 = 0.0;
    for(auto m = 0; m < 100000; ++m) {
      auto ntry    = 0;
      const auto d = qed::iso_dir(rng);
      const auto P = qed::ann::sample(0.3 * d, (-0.3) * d, rng, ntry);
      expect(close(toolbox::norm(P.k0), std::sqrt(1.09), 1e-12));
      expect(close(toolbox::norm(P.k0 + P.k1), 0.0, 0.0, 1e-12));
      const auto Q = qed::ann::sample(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), rng, ntry);
      expect(close(toolbox::norm(Q.k0), 1.0, 1e-12));
      c2 += Q.k0[2] * Q.k0[2];
    }
    expect(close(c2 / 100000, 1.0 / 3.0, 0.02));
  };

  "photon-photon sampler conserves four-momentum"_test = [] {
    auto rng = runko::Rng { runko::rng_key(7, 0, 0, 0), 0 };
    auto n   = 0;
    while(n < 20000) {
      const auto k1 = std::exp(-1.0 + 4.0 * rng.uniform()) * qed::iso_dir(rng);
      const auto k2 = std::exp(-1.0 + 4.0 * rng.uniform()) * qed::iso_dir(rng);
      if(qed::gg::rate(k1, k2) <= 0.0) { continue; }
      ++n;
      auto ntry    = 0;
      const auto P = qed::gg::sample(k1, k2, rng, ntry);
      expect(close(gam(P.p1) + gam(P.p2), toolbox::norm(k1) + toolbox::norm(k2), 1e-9));
      const auto p0 = k1 + k2, p1 = P.p1 + P.p2;
      for(auto i = 0uz; i < 3uz; ++i) { expect(close(p1[i], p0[i], 1e-9, 1e-12)); }
    }
  };

  "bin bounds dominate the rate"_test = [] {
    auto rng = runko::Rng { runko::rng_key(8, 0, 0, 0), 0 };
    for(auto b = 0; b < 2000; ++b) {
      const auto xlo = std::exp(-8.0 + 10.0 * rng.uniform());
      const auto xhi = xlo * std::exp(rng.uniform() * 2.0);
      const auto glo = 1.0 + 30.0 * rng.uniform() * rng.uniform();
      const auto ghi = glo * std::exp(rng.uniform());
      const auto x2lo = std::exp(-8.0 + 10.0 * rng.uniform());
      const auto x2hi = x2lo * std::exp(rng.uniform() * 2.0);
      const auto g2lo = 1.0 + 30.0 * rng.uniform() * rng.uniform();
      const auto g2hi = g2lo * std::exp(rng.uniform());
      const auto bc = qed::comp::bound(xlo, xhi, glo, ghi);
      const auto ba = qed::ann::bound(glo, ghi, g2lo, g2hi);
      const auto bg = qed::gg::bound(xlo, xhi, x2lo, x2hi);
      for(auto m = 0; m < 50; ++m) {
        const auto x  = xlo * std::pow(xhi / xlo, rng.uniform());
        const auto g  = glo * std::pow(ghi / glo, rng.uniform());
        const auto x2 = x2lo * std::pow(x2hi / x2lo, rng.uniform());
        const auto g2 = g2lo * std::pow(g2hi / g2lo, rng.uniform());
        const auto k  = x * qed::iso_dir(rng);
        const auto u  = std::sqrt(g * g - 1.0) * qed::iso_dir(rng);
        const auto k2 = x2 * qed::iso_dir(rng);
        const auto u2 = std::sqrt(g2 * g2 - 1.0) * qed::iso_dir(rng);
        expect(qed::comp::rate(k, u) <= bc * (1.0 + 1e-12)) << "compton";
        expect(qed::ann::rate(u, u2) <= ba * (1.0 + 1e-12)) << "ann";
        expect(qed::gg::rate(k, k2) <= bg * (1.0 + 1e-12)) << "gg";
      }
    }
    // an empty bin (emin = inf, emax = 0) bounds to exactly zero
    const auto inf = std::numeric_limits<double>::infinity();
    expect(qed::comp::bound(inf, 0.0, 1.0, 2.0) == 0.0);
    expect(qed::ann::bound(1.0, 2.0, inf, 0.0) == 0.0);
    expect(qed::gg::bound(inf, 0.0, inf, 0.0) == 0.0);
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{
  return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv }));
}
