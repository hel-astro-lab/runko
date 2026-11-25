// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The channels of the sampler. A two-body channel is a record type the kernels are
// templates over: its two species, the keys of its two sides, what an event does to
// the partner, its rate function and the bound of it over a pair of energy bins.

#include "runko/qed/field.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/qed/species.h"
#include "runko/tools/math.h"

#include <cstddef>
#include <stdexcept>

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

/// The keys one side of a channel spans in every (cell, bin) group: nranges
/// consecutive keys from stride * group + first.
struct [[nodiscard]] Side {
  int stride, first, nranges;

  /// The one key of a species without sub-ranges (the leptons).
  [[nodiscard]] static constexpr Side whole() { return { 1, 0, 1 }; }
  /// One photon sub-range.
  [[nodiscard]] static constexpr Side bit(const int b) { return { 2, b, 1 }; }
  /// Both photon sub-ranges as one range.
  [[nodiscard]] static constexpr Side both_bits() { return { 2, 0, 2 }; }
};

/// One side of a two-body channel as the kernels read it: the one or two sub-ranges
/// of every (cell, bin) seen as one range, and the bounds over them.
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
// One type per two-body channel; the kernels of the pass are templates over them.
//   t1, t2            incident and partner species
//   side1, side2      the keys of each side
//   c_l               pairs counted per slot pair (2 where a species is split in two)
//   consumes_partner  an event removes the partner: lanes may not share one
//   merges            a heavy incident and a light partner merge events (dgmax, g_min)
//   makes_pairs       the event creates lepton pairs at the lepton weight
//   rate, bound       the rate function sigma-hat F / 2 and its bound over two bins

/// Compton scattering of a photon on an electron or a positron; bounded by sigma_KN
/// at the top of the bins.
struct ChCompton {
  Species t2;
  static constexpr Species t1            = Species::ph;
  static constexpr Side side1            = Side::both_bits();
  static constexpr Side side2            = Side::whole();
  static constexpr double c_l            = 1.0;
  static constexpr bool consumes_partner = false;
  static constexpr bool merges           = true;
  static constexpr bool makes_pairs      = false;
  static constexpr double g_min          = 1.0;

  constexpr double rate(const Vec3& k, const Vec3& u, const FieldView&) const
  { return comp::rate(k, u); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return comp::bound(lo1, hi1, lo2, hi2); }
  /// Most energy one scattering of a photon x gives the lepton g.
  constexpr double dgmax(double x, double g) const { return comp::dgmax(x, g); }
};

/// Annihilation of an electron and a positron into two photons; bounded by the
/// annihilation bracket.
struct ChAnn {
  static constexpr Species t1            = Species::em;
  static constexpr Species t2            = Species::ep;
  static constexpr Side side1            = Side::whole();
  static constexpr Side side2            = Side::whole();
  static constexpr double c_l            = 1.0;
  static constexpr bool consumes_partner = true;
  static constexpr bool merges           = false;
  static constexpr bool makes_pairs      = false;

  constexpr double rate(const Vec3& um, const Vec3& up, const FieldView&) const
  { return ann::rate(um, up); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return ann::bound(lo1, hi1, lo2, hi2); }
};

/// Pair production by two photons, one of each sub-range; bounded by sigma_gg.
struct ChGG {
  static constexpr Species t1            = Species::ph;
  static constexpr Species t2            = Species::ph;
  static constexpr Side side1            = Side::bit(0);
  static constexpr Side side2            = Side::bit(1);
  static constexpr double c_l            = 2.0;
  static constexpr bool consumes_partner = true;
  static constexpr bool merges           = false;
  static constexpr bool makes_pairs      = true;

  constexpr double rate(const Vec3& k1, const Vec3& k2, const FieldView&) const
  { return gg::rate(k1, k2); }
  constexpr double bound(double lo1, double hi1, double lo2, double hi2) const
  { return gg::bound(lo1, hi1, lo2, hi2); }
};

/// Synchrotron self-absorption of a photon by an electron or a positron, at the tile
/// field across the lepton's pitch angle; bounded by a scan of the rate over the
/// lepton bin at the photon bin's floor.
struct ChSSA {
  Species t2;
  double b_tile;  // tile field, an upper bound of every lepton's field
  double x_min;   // floor of the synchrotron band

  static constexpr Species t1            = Species::ph;
  static constexpr Side side1            = Side::both_bits();
  static constexpr Side side2            = Side::whole();
  static constexpr double c_l            = 1.0;
  static constexpr bool consumes_partner = false;
  static constexpr bool merges           = true;
  static constexpr bool makes_pairs      = false;
  static constexpr double g_min          = ssa::GMIN;  // no absorption below it

  constexpr double rate(const Vec3& k, const Vec3& u, const FieldView& fv) const
  {
    return ssa::rate(k, u, fv.b_ssa(u / sstd::max(toolbox::norm(u), 1e-300)), x_min);
  }
  constexpr double bound(double lo1, double /*hi1*/, double lo2, double hi2) const
  { return ssa::bound(b_tile, x_min, lo1, lo2, hi2); }
  /// Most energy one absorption of a photon x gives the lepton.
  constexpr double dgmax(double x, double) const { return ssa::dgmax(x, 0.0); }
};

/// Calls f with the record of the two-body channel `which`.
template<typename F>
constexpr auto
  with_record(const Channel which, const double b_tile, const double x_min, F&& f)
{
  switch(which) {
    case Channel::compton_m: return f(ChCompton { Species::em });
    case Channel::compton_p: return f(ChCompton { Species::ep });
    case Channel::pairann: return f(ChAnn {});
    case Channel::photann: return f(ChGG {});
    case Channel::ssa_m: return f(ChSSA { Species::em, b_tile, x_min });
    case Channel::ssa_p: return f(ChSSA { Species::ep, b_tile, x_min });
    default: throw std::runtime_error { "not a two-body channel" };
  }
}

}  // namespace qed
