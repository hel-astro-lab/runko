// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/mdgrid_common.h"
#include "runko/pic/particle.h"
#include "runko/tools/math.h"
#include "thrust/device_vector.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <cstdint>
#include <vector>

namespace qed {

/// Energy in the container's precision: |k| for photons, gamma for leptons.
///
/// Layout, kernels and placement all bin by it, so a slot has one bin everywhere.
[[nodiscard]]
constexpr float
  energy_of(const float mass, const float vx, const float vy, const float vz)
{
  return sstd::sqrt(mass * mass + vx * vx + vy * vy + vz * vz);
}

/// Binning and slack rule of an epoch.
struct LayoutParams {
  /// Energy bins per cell and species.
  int nbin;
  /// Reaction sub-cells per tile edge.
  int ncell_edge;
  /// Slack slots per key, as a fraction of its live LPs.
  double slack;
  /// Floor on the slack slots of a key.
  int smin;
};

/// Device-side tables of a layout, as the kernels read them.
struct LayoutView {
  using I = runko::index_t;
  /// Key of every slot, -1 in the tail.
  std::int32_t* skey;
  /// Per key: first slot.
  I* b;
  /// Per key: slots.
  I* cap;
  /// Per key: live + appended.
  I* fill;
  /// Per key: products appended in the epoch, fitted or not.
  I* prod;
  /// Per key: largest weight held.
  float* wmax;
  /// Per key: lowest energy held.
  float* emin;
  /// Per key: highest energy held.
  float* emax;
  /// Energy bins per cell.
  int nbin;
  /// Reaction sub-cells per tile edge.
  int ncell_edge;
  int ncell;
  int nkeys;
  /// Photons: the key carries the sub-range bit.
  bool keyed;
  float mass;
  /// Log of the lowest kinetic energy of the epoch.
  double llo;
  /// Log width of an energy bin.
  double dl;
  /// Tile extent in global coordinates.
  std::array<float, 3> mins;
  std::array<float, 3> maxs;

  /// Reaction cell of a position.
  [[nodiscard]] constexpr int cell_of(const float x, const float y, const float z) const
  {
    const auto at = [&](const float v, const std::size_t i) {
      const auto f =
        (v - mins[i]) / (maxs[i] - mins[i]) * static_cast<float>(ncell_edge);
      return sstd::clamp(static_cast<int>(sstd::floor(f)), 0, ncell_edge - 1);
    };
    return (at(x, 0uz) * ncell_edge + at(y, 1uz)) * ncell_edge + at(z, 2uz);
  }

  /// Energy bin of energy e, log-spaced in the kinetic energy of the epoch.
  [[nodiscard]] constexpr int bin_of(const double e) const
  {
    if(nbin == 1) { return 0; }
    const auto v = sstd::max(e - static_cast<double>(mass), 1e-30);
    return sstd::clamp(
      static_cast<int>(sstd::floor((sstd::log(v) - llo) / dl)),
      0,
      nbin - 1);
  }

  /// Key of a (cell, bin) pair and photon bit.
  [[nodiscard]] constexpr int
    key_of(const int cell, const int ebin, const int bit) const
  {
    const auto k = cell * nbin + ebin;
    return keyed ? 2 * k + bit : k;
  }

  /// Energy bin of a key.
  [[nodiscard]] constexpr int bin_of_key(const int k) const
  {
    return (keyed ? k / 2 : k) % nbin;
  }

  /// Widens the bounds of key k over a member of energy e and weight w.
  constexpr void widen(const int k, const float e, const float w) const
  {
    sstd::atomic_max(&wmax[k], w);
    sstd::atomic_max(&emax[k], e);
    sstd::atomic_min(&emin[k], e);
  }
};

/// Weight of each LP, payload component 0, of a species with payload (the photons);
/// an empty view for one without (the leptons).
template<typename Container>
[[nodiscard]] auto
  weight_mds(Container& con)
{
  using M = decltype(con.payload_mds(0));
  return con.n_payload() > 0uz ? con.payload_mds(0) : M {};
}

/// Padded, binned layout of one species.
///
/// The slots [b_k, b_k + cap_k) of key k hold its live LPs, then slack for products;
/// slots past the last key are the tail. The key is the reaction cell, the energy bin
/// and, for photons, a sub-range bit that makes photon-photon pair production a
/// cross-range channel. Slack slots are dead ids; a product written there becomes
/// alive. The container's append fills the dead slots behind its last alive one, so
/// an arrival may sit in a key's slack: a product claims a slack slot only while its
/// id is still dead.
///
/// The container belongs to the tile and is passed to every call, so a tile may be
/// copied or moved without a stale pointer.
class [[nodiscard]] Layout {
  /// Photons carry the sub-range bit.
  bool keyed_;
  /// Rest mass: 0 photons, 1 leptons.
  float mass_;
  /// LP weight without payload: leptons.
  double w0_;
  /// Parameters of the current epoch.
  LayoutParams par_ {};
  std::array<float, 3> mins_ {}, maxs_ {};
  /// Log bin edges: lowest edge and bin width.
  double llo_ { 0.0 }, dl_ { 1.0 };
  int nkeys_ { 0 };
  bool laid_ { false };
  /// The container's reorder count at the build.
  std::uint64_t reorders_ { 0 };

  // The tables of LayoutView.
  thrust::device_vector<std::int32_t> skey_;
  thrust::device_vector<runko::index_t> b_, cap_, fill_, prod_;
  thrust::device_vector<float> wmax_, emin_, emax_;

  /// Layout scratch: per old slot its key and the sort permutation, per key + 1 the
  /// key starts, per new slot its source; the high-water mark is kept between epochs.
  thrust::device_vector<runko::index_t> keys_, ord_, cnt_, src_;

public:
  /// The weight of every LP of a species without payload is w0; a species with one
  /// carries the weight of each LP as payload component 0.
  Layout(bool keyed, double mass, double w0);

  [[nodiscard]] double w0() const { return w0_; }

  /// Whether the tables describe con: built, not reordered since (a sort drops the
  /// layout), and at the size they were built or resized for.
  [[nodiscard]] bool laid(const pic::ParticleContainer& con) const
  {
    return laid_ and con.reorder_count() == reorders_ and skey_.size() == con.size();
  }
  void invalidate() { laid_ = false; }

  /// Groups the live LPs by key with slack behind each group (one reorder of the
  /// container) and builds the tables.
  void build(
    pic::ParticleContainer& con,
    const LayoutParams& par,
    std::array<float, 3> mins,
    std::array<float, 3> maxs,
    std::uint64_t rng_key);

  /// The tables as the kernels read them: raw device pointers, valid while laid and
  /// until the next resize.
  [[nodiscard]] LayoutView view();

  /// Follows the container's size since the last build or resize.
  ///
  /// Slots it grew by (products placed, arrivals appended) are tail (-1). Slots it
  /// shrank by were dead slack (an append reuses the dead slots behind the last alive
  /// one and drops the rest); they are appended again so that every key keeps its
  /// slots.
  void resized(pic::ParticleContainer& con);

  /// Widens the bounds of every key by the live LPs it holds, on the caller's queue;
  /// the caller waits.
  void touch_all(pic::ParticleContainer& con, const tyvi::mdgrid_work& w);

  /// Live LPs and their summed weight, energy and momentum.
  struct Totals {
    double N, W, E;
    std::array<double, 3> P;
  };
  [[nodiscard]] Totals totals(const pic::ParticleContainer& con) const;

  /// Host copy of the key of every slot, for IO and tests.
  [[nodiscard]] std::vector<std::int32_t>
    get_slot_keys(const pic::ParticleContainer& con);

  /// Host copy of the per-key tables and the bin edges, for IO and tests.
  struct Tables {
    std::vector<runko::index_t> b, cap, fill, prod;
    std::vector<float> wmax, emin, emax;
    double llo, dl;
  };
  [[nodiscard]] Tables get_tables();
};

}  // namespace qed
