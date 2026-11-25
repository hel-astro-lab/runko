// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Padded, binned layout of one species: the slots [b_k, b_k + cap_k) of key k hold its
// live LPs, then slack for products; slots past the last key are the tail. The key is
// the reaction cell, the energy bin and, for photons, a sub-range bit that makes
// photon-photon pair production a cross-range channel. Slack slots are dead ids; a
// product written there becomes alive. The container's append fills the dead slots
// behind its last alive one, so an arrival may sit in a key's slack: a product claims a
// slack slot only while its id is still dead.

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
/// Layout, kernels and placement all bin by it, so a slot has one bin everywhere.
[[nodiscard]]
constexpr float
  energy_of(const float mass, const float vx, const float vy, const float vz)
{ return sstd::sqrt(mass * mass + vx * vx + vy * vy + vz * vz); }

/// Binning and slack rule of an epoch.
struct [[nodiscard]] LayoutParams {
  int nbin;        // energy bins per cell and species
  int ncell_edge;  // reaction sub-cells per tile edge
  double slack;    // slack slots per key, as a fraction of its live LPs
  int smin;        // floor on the slack slots of a key
};

/// Device-side tables of a layout, as the kernels read them.
struct [[nodiscard]] LayoutView {
  using I = runko::index_t;
  std::int32_t* skey;  // key of every slot, -1 in the tail
  I* b;                // per key: first slot
  I* cap;              // per key: slots
  I* fill;             // per key: live + appended
  I* prod;             // per key: products appended in the epoch, fitted or not
  float* wmax;         // per key: largest weight held
  float* emin;         // per key: lowest energy held
  float* emax;         // per key: highest energy held
  int nbin;            // energy bins per cell
  int ncell_edge;      // reaction sub-cells per tile edge
  int ncell;           // reaction sub-cells of the tile
  int nkeys;
  bool keyed;                 // photons: the key carries the sub-range bit
  float mass;                 // rest mass of the species
  double llo;                 // log of the lowest kinetic energy of the epoch
  double dl;                  // log width of an energy bin
  std::array<float, 3> mins;  // tile extent in global coordinates
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
  { return (keyed ? k / 2 : k) % nbin; }
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

/// Layout of one species. The container belongs to the tile and is passed to every
/// call, so a tile may be copied or moved without a stale pointer.
class [[nodiscard]] Pop {
  bool keyed_;                              // photons carry the sub-range bit
  float mass_;                              // rest mass: 0 photons, 1 leptons
  double w0_;                               // LP weight without payload: leptons
  LayoutParams par_ {};                     // of the current epoch
  std::array<float, 3> mins_ {}, maxs_ {};  // tile extent
  double llo_ { 0.0 }, dl_ { 1.0 };         // log bin edges: lowest edge and bin width
  int nkeys_ { 0 };
  bool laid_ { false };           // tables built
  std::uint64_t reorders_ { 0 };  // the container's reorder count at the build

  // the tables of LayoutView
  thrust::device_vector<std::int32_t> skey_;
  thrust::device_vector<runko::index_t> b_, cap_, fill_, prod_;
  thrust::device_vector<float> wmax_, emin_, emax_;

  // layout scratch: per old slot its key and the sort permutation, per key + 1 the key
  // starts, per new slot its source; the high-water mark is kept between epochs
  thrust::device_vector<runko::index_t> keys_, ord_, cnt_, src_;


public:
  /// w0 is the weight of every LP of a species without payload; a species with one
  /// carries the weight of each LP as payload component 0.
  Pop(bool keyed, double mass, double w0);

  [[nodiscard]] double w0() const { return w0_; }

  /// Whether the tables describe con: built, not reordered since (a sort drops the
  /// layout), and at the size they were built or resized for.
  [[nodiscard]] bool laid(const pic::ParticleContainer& con) const
  { return laid_ and con.reorder_count() == reorders_ and skey_.size() == con.size(); }
  void invalidate() { laid_ = false; }

  /// Groups the live LPs by key with slack behind each group (one reorder of the
  /// container) and builds the tables.
  void layout(
    pic::ParticleContainer& con,
    const LayoutParams& par,
    std::array<float, 3> mins,
    std::array<float, 3> maxs,
    std::uint64_t rng_key);

  /// The tables as the kernels read them: raw device pointers, valid while laid and
  /// until the next resize.
  [[nodiscard]] LayoutView view();

  /// Follows the container's size since the last build or resize. Slots it grew by
  /// (products placed, arrivals appended) are tail (-1). Slots it shrank by were dead
  /// slack (an append reuses the dead slots behind the last alive one and drops the
  /// rest); they are appended again so that every key keeps its slots.
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

  // ---- host copies for IO and tests

  [[nodiscard]] std::vector<std::int32_t>
    get_slot_keys(const pic::ParticleContainer& con);

  struct Tables {
    std::vector<runko::index_t> b, cap, fill, prod;  // per key, as on the device
    std::vector<float> wmax, emin, emax;             // per key
    double llo, dl;                                  // log bin edges
  };
  [[nodiscard]] Tables get_tables();
};

}  // namespace qed
