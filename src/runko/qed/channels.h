// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The channels of the sampler. A two-body channel is a record type the kernels are
// templates over: its two species, how an event treats the partner, the slot ranges of
// its two sides, its rate function and the bound of it over a pair of energy bins.

#include "runko/qed/field.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/qed/species.h"
#include "runko/tools/math.h"

#include <cstddef>

namespace qed {

/// Every channel; the order numbers the random streams of the kernels.
enum class Channel : int {
  compton_m = 0,
  compton_p,
  pairann,
  photann,
  ssa_m,
  ssa_p,
  syn_m,
  syn_p,
  bw,
  count
};
inline constexpr auto N_CHANNELS = static_cast<std::size_t>(Channel::count);

[[nodiscard]] const char* channel_name(Channel c);
[[nodiscard]] bool is_twobody(Channel c);

/// What an event does to the partner: kicks it, consumes it, or feeds it energy.
enum class Kind { scatter, convert, absorb };

/// The slot ranges a side spans: the cell's range, one photon bit, or both bits.
enum class Sel { cell, bit0, bit1, both };

/// One side of a two-body channel: the one or two sub-ranges of every (cell, bin) seen
/// as one range.
struct [[nodiscard]] SideView {
  using I = runko::index_t;
  const I* base0;     // first slot of sub-range 0, per (cell, bin)
  const I* cnt0;      // its slots (live + appended)
  const I* base1;     // first slot of sub-range 1 (second photon bit), 0 if absent
  const I* cnt1;      // its slots, 0 if absent
  const float* wmax;  // largest weight over both sub-ranges
  const float* emin;  // energy bounds over both sub-ranges
  const float* emax;

  [[nodiscard]] constexpr I N(const int cb) const { return cnt0[cb] + cnt1[cb]; }

  /// Slot of the virtual index v in [0, N) of group cb.
  [[nodiscard]] constexpr I pick(const int cb, const I v) const
  {
    const auto c0 = cnt0[cb];
    return v < c0 ? base0[cb] + v : base1[cb] + v - c0;
  }
};

// ---- the two-body records
// One type per two-body channel; the sampler kernels are templates over them.

// Compton: photon on e- or e+, scatters; bounded by sigma_KN at the top of the bins.
struct ChCompton {
  Species t2;
  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::scatter;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  static constexpr bool pairs   = false;
  static constexpr double g_min = 1.0;

  constexpr double rate(const Vec3& k, const Vec3& u) const { return comp::rate(k, u); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return comp::bound(lo1, hi1, lo2, hi2); }
  constexpr double dgmax(double x, double g) const { return comp::dgmax(x, g); }
};

// Pair annihilation: e- on e+, both convert into photons; bounded by the annihilation
// bracket.
struct ChAnn {
  static constexpr Species t1   = Species::em;
  static constexpr Species t2   = Species::ep;
  static constexpr Kind kind    = Kind::convert;
  static constexpr Sel sel1     = Sel::cell;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = false;
  static constexpr bool pairs   = false;

  constexpr double rate(const Vec3& um, const Vec3& up) const
  { return ann::rate(um, up); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return ann::bound(lo1, hi1, lo2, hi2); }
};

// Photon-photon pair production: photon bit 0 on photon bit 1, converts; bounded by
// sigma_gg.
struct ChGG {
  static constexpr Species t1   = Species::ph;
  static constexpr Species t2   = Species::ph;
  static constexpr Kind kind    = Kind::convert;
  static constexpr Sel sel1     = Sel::bit0;
  static constexpr Sel sel2     = Sel::bit1;
  static constexpr double c_l   = 2.0;  // same species split by the photon bit
  static constexpr bool has_dgm = false;
  static constexpr bool pairs   = true;

  constexpr double rate(const Vec3& k1, const Vec3& k2) const
  { return gg::rate(k1, k2); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return gg::bound(lo1, hi1, lo2, hi2); }
};

// Synchrotron self-absorption: photon on e- or e+, the lepton absorbs it; bounded by a
// scan of the rate over the lepton bin at the photon bin's floor.
struct ChSSA {
  Species t2;
  double b_tile;  // tile field, an upper bound of every lepton's b_ssa
  double x_min;   // floor of the synchrotron band

  static constexpr Species t1   = Species::ph;
  static constexpr Kind kind    = Kind::absorb;
  static constexpr Sel sel1     = Sel::both;
  static constexpr Sel sel2     = Sel::cell;
  static constexpr double c_l   = 1.0;
  static constexpr bool has_dgm = true;
  static constexpr bool pairs   = false;
  static constexpr double g_min = ssa::GMIN;  // no absorption below the threshold

  /// Rate at the lepton's field b, set by its pitch angle in the tile field.
  constexpr double rate(const Vec3& k, const Vec3& u, const double b) const
  { return ssa::rate(k, u, b, x_min); }

  constexpr double bound(double lo1, double /*hi1*/, double lo2, double hi2) const
  { return ssa::bound(b_tile, x_min, lo1, lo2, hi2); }

  constexpr double dgmax(double x, double) const { return ssa::dgmax(x, 0.0); }
};

/// Rate of the pair (p1, p2).
template<typename Ch>
[[nodiscard]]
constexpr double
  pair_rate(const Ch& ch, const FieldView& fv, const Vec3& p1, const Vec3& p2)
{
  if constexpr(Ch::kind == Kind::absorb) {
    return ch.rate(p1, p2, fv.b_ssa(p2 / sstd::max(toolbox::norm(p2), 1e-300)));
  } else {
    return ch.rate(p1, p2);
  }
}

}  // namespace qed
