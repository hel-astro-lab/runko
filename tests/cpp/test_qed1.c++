// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/tools/random.h"

#include <boost/ut.hpp>
#include <cmath>

namespace {

using namespace boost::ut;
using qed::Vec3;

constexpr bool
  close(const double a, const double b, const double rtol, const double atol = 0.0)
{ return std::abs(a - b) <= rtol * std::abs(b) + atol; }

[[maybe_unused]]
const suite<"one-body and self-absorption physics"> _ = [] {
  "synchrotron pins"_test = [] {
    expect(close(qed::syn::T(1e-6), 1.4433737800784523, 1e-11));
    expect(close(qed::syn::T(1e-3), 1.445154513321066, 1e-11));
    expect(close(qed::syn::T(0.1), 1.3474022710638174, 1e-11));
    expect(close(qed::syn::T(1.0), 1.0291979856368023, 1e-11));
    expect(close(qed::syn::T(10.0), 0.6141058186542885, 1e-11));
    expect(close(qed::syn::T(1e3), 0.1433893703388872, 1e-11));

    expect(close(qed::syn::g(1e-3), 0.99459431527274, 1e-11));
    expect(close(qed::syn::g(0.1), 0.6627611620800354, 1e-11));
    expect(close(qed::syn::g(1.0), 0.18110213894202243, 1e-11));
    expect(close(qed::syn::g(10.0), 0.01849384790421692, 1e-11));

    expect(close(qed::syn::P(1e-3, 1e-4), 0.474697273534346, 1e-10));
    expect(close(qed::syn::P(1e-3, 1e-2), 0.9987942542651651, 1e-10));
    expect(close(qed::syn::P(0.1, 0.05), 0.789712828904259, 1e-10));
    expect(close(qed::syn::P(1.0, 0.5), 0.9505199654008519, 1e-10));
    expect(close(qed::syn::P(10.0, 0.9), 0.9794974779417613, 1e-10));

    expect(close(qed::syn::r_of(0.1, 1e-3), 7.591896296745883e-07, 1e-8));
    expect(close(qed::syn::r_of(0.5, 1e-3), 0.00011888364393755561, 1e-8));
    expect(close(qed::syn::r_of(0.9, 1e-3), 0.0013536481304822032, 1e-8));
    expect(close(qed::syn::r_of(0.3, 0.1), 0.0018254245223500013, 1e-8));
    expect(close(qed::syn::r_of(0.7, 1.0), 0.1302365994192557, 1e-8));
    expect(close(qed::syn::r_of(0.999, 0.05), 0.28903472579233747, 1e-8));
    expect(close(qed::syn::r_of(0.999, 5.0), 0.9858025878137909, 1e-8));

    expect(close(qed::syn::rmean(1e-3, 0.0), 0.0004915576082575523, 1e-9));
    expect(close(qed::syn::rmean(0.1, 0.0), 0.03293233949215805, 1e-9));
    expect(close(qed::syn::rmean(1.0, 0.0), 0.1190278525361551, 1e-9));
    expect(close(qed::syn::rmean(1e-3, 1e-4), 0.0009150840279317343, 1e-9));
    expect(close(qed::syn::rmean(0.1, 0.01), 0.0648880516895646, 1e-9));

    const auto u = Vec3(3.0, 1.0, -2.0);
    expect(close(qed::syn::rate(u, 1e-3, 0.0, false), 0.0014452134186639385, 1e-10));
    expect(close(qed::syn::rate(u, 1e-3, 1e-3, false), 0.0008376097799749066, 1e-10));
    expect(qed::syn::rate(Vec3(2.0, 0.0, 0.0), 1e-3, 0.0, false) == 0.0);  // gamma < 3
  };

  "the closed-form mean above a floor follows the moment of the cumulative"_test = [] {
    // reference: Simpson's rule in ln r on the cumulative's tail from a few decades
    // below the spectrum's scale; the sliver below the first node in closed form
    const auto ref = [](const double chi, const double rmin) {
      const auto s    = qed::syn::shape(chi);
      const auto P0   = qed::syn::P(s, rmin);
      auto r_lo       = std::max(rmin, 1e-7 * std::min(chi, 1.0));
      r_lo            = std::min(r_lo, 1.0 - 1e-12);
      constexpr int n = 4096;
      const auto lr0  = std::log(r_lo);
      auto sum        = 0.0;
      for(auto i = 0; i <= n; ++i) {
        const auto r = std::exp(lr0 * (1.0 - static_cast<double>(i) / n));
        const auto w = (i == 0 or i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        sum += w * (1.0 - qed::syn::P(s, r)) * r;
      }
      const auto h = -lr0 / n;
      return rmin + (h / 3.0 * sum + (r_lo - rmin) * (1.0 - qed::syn::P(s, r_lo))) /
                      (1.0 - P0);
    };

    // the fitted range, every floor that leaves the channel open (1 - P0 >= FMIN)
    for(const auto chi: { 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0 }) {
      const auto s   = qed::syn::shape(chi);
      const auto tol = chi <= 1e-2 ? 3e-3 : (chi <= 1.0 ? 5e-3 : 8e-3);
      for(auto z = 0.0; z <= 9.0; z += 0.5) {
        // the floor of tail exponent z
        const auto rmin = std::pow(1.0 - std::exp(-z / (4.0 * s.n)), 3.0);
        if(qed::syn::frac(chi, rmin) <= 0.0) { continue; }
        expect(close(qed::syn::rmean(s, rmin), ref(chi, rmin), tol))
          << "chi " << chi << " z " << z;
      }
    }
  };

  "synchrotron inversion round-trips and the sampler conserves energy"_test = [] {
    auto rng = runko::Rng { runko::rng_key(11, 0, 0, 0), 0 };

    // round-trip to round-off except where the quartic degenerates (e.g. chi ~ 2.6, as
    // in the prototype): those count as outliers
    auto bad = 0;
    for(auto n = 0; n < 20000; ++n) {
      const auto chi = std::exp(-12.0 + 16.0 * rng.uniform());
      const auto z   = qed::syn::ZCLIP + (1.0 - 2.0 * qed::syn::ZCLIP) * rng.uniform();
      const auto r   = qed::syn::r_of(z, chi);
      expect(r >= 0.0 and r <= 1.0);
      if(r > 0.0 and r < 1.0) { bad += std::abs(qed::syn::P(chi, r) - z) > 1e-4; }
    }
    expect(bad < 40) << bad;

    auto worst = 1.0, below = 0.0, open = 0.0;
    for(auto n = 0; n < 20000; ++n) {
      const auto gam = 3.0 + 100.0 * rng.uniform();
      const auto u   = std::sqrt(gam * gam - 1.0) * qed::iso_dir(rng);
      const auto b   = std::exp(-8.0 + 6.0 * rng.uniform());
      const auto xm  = 0.5 * rng.uniform() * gam * rng.uniform();
      const auto P   = qed::syn::sample(u, b, rng, false, xm);
      const auto x   = toolbox::norm(P.k);
      const auto g1  = std::sqrt(1.0 + toolbox::dot(P.u1, P.u1));
      expect(x >= 0.0);
      if(x > 0.0) {
        expect(close(toolbox::dot(P.k, u) / (x * toolbox::norm(u)), 1.0, 1e-9))
          << "collinear";
      }
      expect(close(g1 + x, gam, 1e-9) or g1 == 1.0) << "energy";
      // where the channel is open, the floor holds up to the round-off of the inverted
      // Pade cumulative
      if(qed::syn::rate(u, b, xm, false) > 0.0) {
        open += 1.0;
        worst = std::min(worst, x / xm);
        below += (x < xm);
      }
    }

    expect(open > 5000.0);
    expect(below / open < 0.02) << below / open;
    expect(worst > 0.9) << worst;
  };

  "breit-wheeler pins"_test = [] {
    expect(close(qed::bw::T(0.1), 5.922685739450321e-13, 1e-10));
    expect(close(qed::bw::T(0.5), 0.0010315347708528574, 1e-11));
    expect(close(qed::bw::T(1.0), 0.014130758374273423, 1e-11));
    expect(close(qed::bw::T(10.0), 0.10848430195842701, 1e-11));
    expect(close(qed::bw::T(1e3), 0.03729284600777369, 1e-11));

    expect(close(qed::bw::P(0.5, 0.6), 0.7578413698673724, 1e-10));
    expect(close(qed::bw::P(1.0, 0.9), 0.9987354700762969, 1e-10));
    expect(close(qed::bw::P(10.0, 0.75), 0.7446145277556716, 1e-10));
    // clamped to chi = 0.5
    expect(close(qed::bw::P(0.3, 0.6), 0.7578413698673724, 1e-10));
    expect(close(qed::bw::P(1.0, 0.2), 0.029197882807714626, 1e-10));

    expect(close(qed::bw::r_of(0.6, 0.5), 0.5362842612233096, 1e-8));
    expect(close(qed::bw::r_of(0.95, 1.0), 0.772002212294584, 1e-8));
    expect(close(qed::bw::r_of(0.3, 10.0), 0.2919230665454289, 1e-8));
    expect(close(qed::bw::r_of(0.5, 3.0), 0.5, 1e-8));
    expect(close(qed::bw::r_of(0.999, 100.0), 0.998501706239715, 1e-8));

    expect(close(qed::bw::rate(Vec3(5.0, 0.0, 0.0), 0.3), 0.009896635737587917, 1e-10));
    expect(qed::bw::rate(Vec3(1.5, 0.0, 0.0), 0.3) == 0.0);
  };

  "breit-wheeler symmetry, round-trip and energy conservation"_test = [] {
    auto rng = runko::Rng { runko::rng_key(12, 0, 0, 0), 0 };

    for(auto n = 0; n < 20000; ++n) {
      const auto chi = std::exp(-2.0 + 12.0 * rng.uniform());
      const auto z   = qed::bw::ZCLIP + (1.0 - 2.0 * qed::bw::ZCLIP) * rng.uniform();
      const auto r   = qed::bw::r_of(z, chi);
      expect(close(r + qed::bw::r_of(1.0 - z, chi), 1.0, 0.0, 1e-12)) << "reflection";
      expect(close(qed::bw::P(chi, r), z, 0.0, 2e-6)) << chi << " " << z;
    }

    auto rsum = 0.0;
    for(auto n = 0; n < 20000; ++n) {
      const auto x  = 2.0 + 100.0 * rng.uniform();
      const auto k  = x * qed::iso_dir(rng);
      const auto b  = std::exp(-4.0 + 4.0 * rng.uniform());
      const auto P  = qed::bw::sample(k, b, rng);
      const auto gm = std::sqrt(1.0 + toolbox::dot(P.um, P.um));
      const auto gp = std::sqrt(1.0 + toolbox::dot(P.up, P.up));
      expect(close(gm + gp, x, 1e-9)) << "energy";
      rsum += (gm - 1.0) / (x - 2.0);
    }

    expect(close(rsum / 20000, 0.5, 0.01)) << "charge symmetry";
  };

  "bessel fits against scipy"_test = [] {
    expect(close(qed::ssa::int_K53(1e-6), 21493.46861598466, 5e-3));
    expect(close(qed::ssa::int_K53(1e-3), 213.13906509157854, 5e-3));
    expect(close(qed::ssa::int_K53(0.1), 8.181855378460552, 5e-3));
    expect(close(qed::ssa::int_K53(1.0), 0.6514229738248462, 5e-3));
    expect(close(qed::ssa::int_K53(5.0), 0.004249632987541319, 5e-3));
    expect(close(qed::ssa::int_K53(20.0), 5.984363524870444e-10, 5e-3));

    expect(close(qed::ssa::K23(1e-6), 10747.641081108552, 1e-3));
    expect(close(qed::ssa::K23(1e-3), 107.46383549069988, 1e-3));
    expect(close(qed::ssa::K23(0.1), 4.752962677620825, 1e-3));
    expect(close(qed::ssa::K23(1.0), 0.4944750621042118, 1e-3));
    expect(close(qed::ssa::K23(5.0), 0.00384442463449682, 1e-3));
    expect(close(qed::ssa::K23(20.0), 5.803848427192581e-10, 1e-3));

    expect(close(qed::ssa::G_syn(1e-3, 1e-6), 0.0003314650058238534, 5e-3));
    expect(close(qed::ssa::G_syn(1e-3, 1e-4), 2.772040750586381e-35, 2e-2));
    expect(close(qed::ssa::G_syn(0.1, 0.05), 0.00013664151853163116, 5e-3));
    expect(close(qed::ssa::G_syn(1.0, 0.3), 0.28777529652900524, 5e-3));
    expect(close(qed::ssa::G_syn(10.0, 5.0), 1.9911393222639286, 5e-3));
    expect(qed::ssa::G_syn(1.0, 1.0) == 0.0);

    expect(close(qed::b_from_ellB(10.0, 4.111e9), 0.0010000097266020366, 1e-12));
  };

  "self-absorption rate against the tabulated prototype"_test = [] {
    constexpr auto b = 1e-3;
    const auto at    = [&](const double g, const double x) {
      return qed::ssa::rate(
        Vec3(x, 0.0, 0.0),
        Vec3(std::sqrt(g * g - 1.0), 0.0, 0.0),
        b,
        0.0);
    };

    expect(close(at(5.0, 1e-3), 3919.6465542179676, 5e-2)) << at(5.0, 1e-3);
    expect(close(at(20.0, 1e-3), 378.76502927171197, 5e-2)) << at(20.0, 1e-3);
    expect(close(at(20.0, 1e-5), 816049.4545241046, 5e-2)) << at(20.0, 1e-5);
    expect(close(at(100.0, 1e-2), 0.5571555116430739, 5e-2)) << at(100.0, 1e-2);
    expect(close(at(3.0, 1e-4), 451989.0660459265, 5e-2)) << at(3.0, 1e-4);
    expect(at(2.9, 1e-4) == 0.0);

    // the derivative branch joins the difference branch continuously
    for(const auto g: { 5.0, 20.0, 100.0 }) {
      const auto lo = at(g, 0.999 * qed::ssa::XREL * g);
      const auto hi = at(g, 1.001 * qed::ssa::XREL * g);
      expect(close(lo, hi, 2e-2)) << g << " " << lo << " " << hi;
    }

    expect(close(qed::ssa::shmax(b, 1e-3, 3.0, 100.0), 10103.600373061168, 5e-2));
    expect(close(qed::ssa::shmax(b, 1e-2, 3.0, 30.0), 179.75788989184298, 5e-2));
    expect(close(qed::ssa::shmax(b, 1e-5, 10.0, 1000.0), 2735775.4653898813, 5e-2));
  };

  "self-absorption bound dominates the rate; product on shell"_test = [] {
    auto rng = runko::Rng { runko::rng_key(13, 0, 0, 0), 0 };

    for(auto n = 0; n < 300; ++n) {
      const auto b   = std::exp(-9.0 + 5.0 * rng.uniform());
      const auto xlo = std::exp(-12.0 + 8.0 * rng.uniform());
      const auto xhi = xlo * std::exp(3.0 * rng.uniform());
      const auto glo = 3.0 * std::exp(4.0 * rng.uniform());
      const auto ghi = glo * std::exp(1.5 * rng.uniform());
      const auto sh  = qed::ssa::bound(b, 0.0, xlo, glo, ghi);
      for(auto m = 0; m < 200; ++m) {
        const auto x = xlo * std::pow(xhi / xlo, rng.uniform());
        const auto g = glo * std::pow(ghi / glo, rng.uniform());
        const auto k = x * qed::iso_dir(rng);
        const auto u = std::sqrt(g * g - 1.0) * qed::iso_dir(rng);
        expect(qed::ssa::rate(k, u, b, 0.0) <= sh) << b << " " << x << " " << g;
      }
    }

    for(auto n = 0; n < 1000; ++n) {
      const auto g  = 3.0 + 50.0 * rng.uniform();
      const auto u  = std::sqrt(g * g - 1.0) * qed::iso_dir(rng);
      const auto k  = std::exp(-8.0 + 8.0 * rng.uniform()) * qed::iso_dir(rng);
      const auto u1 = qed::ssa::sample(k, u);
      const auto g1 = std::sqrt(1.0 + toolbox::dot(u1, u1));
      expect(close(g1, g + toolbox::norm(k), 1e-12));
    }

    expect(
      qed::ssa::bound(1e-3, 0.0, std::numeric_limits<double>::infinity(), 3.0, 10.0) ==
      0.0);
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{ return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv })); }
