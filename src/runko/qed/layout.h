// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The padded, binned layout of one species: the slots [b_k, b_k + cap_k) of key k hold the live
// LPs of the key first and then slack for products; the slots beyond the last key are the
// tail. The key is the reaction cell, the energy bin and, for photons, the sub-range bit that
// makes photon-photon pair production a cross-range channel. The container itself is the
// layout: slack slots are dead ids, products written into them become alive.

#include "runko/mdgrid_common.h"
#include "runko/pic/particle.h"
#include "runko/tools/math.h"
#include "thrust/device_vector.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <cstdint>
#include <vector>

namespace qed {

/// energy of a momentum in the container's own precision: |k| for photons, gamma for leptons.
/// The one function the layout, the kernels and the placement bin by, so that a slot's bin is
/// the same wherever it is computed.
[[nodiscard]]
constexpr float
  energy_of(const float mass, const float vx, const float vy, const float vz)
{
  return sstd::sqrt(mass * mass + vx * vx + vy * vy + vz * vz);
}

struct [[nodiscard]] LayoutParams {
  int nbin;         // energy bins per cell and species
  int ncell_edge;   // reaction sub-cells per tile edge
  double slack;     // slack slots per key, as a fraction of its live LPs
  int smin;         // and a floor on that count
};

/// The device-side description of a layout: what the kernels read.
struct [[nodiscard]] LayoutView {
  using I = runko::index_t;
  std::int32_t* skey;   // key of every slot, -1 in the tail
  I* b;                 // per key: first slot
  I* cap;               // per key: slots
  I* fill;              // per key: live + appended
  I* prod;              // per key: products appended in the epoch, fitted or not
  float* wmax;          // per key: the bounds the envelopes read
  float* emin;          // per key: lowest energy held
  float* emax;          // per key: highest energy held
  int nbin;             // energy bins per cell
  int ncell_edge;       // reaction sub-cells per tile edge
  int ncell;            // reaction sub-cells of the tile
  int nkeys;            // keys of the layout
  bool keyed;           // photons: the key carries the sub-range bit
  float mass;           // rest mass of the species
  double llo;           // log of the lowest kinetic energy of the epoch
  double dl;            // log width of an energy bin
  std::array<float, 3> mins;  // tile extent, for the reaction cell of a position
  std::array<float, 3> maxs;  // (positions are in the tile's own coordinates)

  /// reaction cell of a position
  [[nodiscard]] constexpr int cell_of(const float x, const float y, const float z) const
  {
    const auto at = [&](const float v, const std::size_t i) {
      const auto f = (v - mins[i]) / (maxs[i] - mins[i]) * static_cast<float>(ncell_edge);
      return sstd::clamp(static_cast<int>(sstd::floor(f)), 0, ncell_edge - 1);
    };
    return (at(x, 0uz) * ncell_edge + at(y, 1uz)) * ncell_edge + at(z, 2uz);
  }

  /// energy bin of an LP of energy e: log-spaced in the kinetic energy of the epoch
  [[nodiscard]] constexpr int bin_of(const double e) const
  {
    if(nbin == 1) { return 0; }
    const auto v = sstd::max(e - static_cast<double>(mass), 1e-30);
    return sstd::clamp(static_cast<int>(sstd::floor((sstd::log(v) - llo) / dl)), 0, nbin - 1);
  }

  /// key of a (cell, bin) pair, and of the photon bit
  [[nodiscard]] constexpr int key_of(const int cell, const int ebin, const int bit) const
  {
    const auto k = cell * nbin + ebin;
    return keyed ? 2 * k + bit : k;
  }

  /// bin of a key
  [[nodiscard]] constexpr int bin_of_key(const int k) const
  {
    return (keyed ? k / 2 : k) % nbin;
  }
};

class [[nodiscard]] Pop {
  // The layout of one species; the container itself belongs to the tile, which passes it to
  // every call, so that a tile may be copied or moved without a pointer going stale.
  // ---- what the species is and what the current epoch fixed
  bool keyed_;                              // photons carry the sub-range bit
  float mass_;                              // rest mass: 0 photons, 1 leptons
  LayoutParams par_ {};                     // bins, cells and slack rule of the epoch
  std::array<float, 3> mins_ {}, maxs_ {};  // tile extent
  double llo_ { 0.0 }, dl_ { 1.0 };         // log bin edges: lowest edge and bin width
  int nkeys_ { 0 };                         // keys of the layout
  bool laid_ { false };                     // tables built; valid while the container keeps their size

  // ---- per slot: key of every slot, -1 in the tail
  thrust::device_vector<std::int32_t> skey_;
  // ---- per key: first slot, slots, live + appended, products appended in the epoch
  thrust::device_vector<runko::index_t> b_, cap_, fill_, prod_;
  // ---- per key: the weight and energy bounds the envelopes read
  thrust::device_vector<float> wmax_, emin_, emax_;
  // ---- scratch of a layout, per slot (key, sort permutation, target slot) and per key + 1
  // (key starts); the high-water mark is kept between epochs
  thrust::device_vector<runko::index_t> keys_, ord_, cnt_, slots_;


public:
  Pop(bool keyed, double mass);

  /// whether the tables describe this container: built, and the container at the size they were
  /// built or grown for (pic's append writes behind the last alive slot and truncates the trailing
  /// slack, so any other size means the tables are stale)
  [[nodiscard]] bool laid(const pic::ParticleContainer& con) const
  {
    return laid_ and skey_.size() == con.size();
  }
  void invalidate() { laid_ = false; }

  /// group the live LPs by key, reserve slack behind each group, build the tables
  void layout(
    pic::ParticleContainer& con,
    const LayoutParams& par,
    std::array<float, 3> mins,
    std::array<float, 3> maxs,
    std::uint64_t rng_key);

  /// the tables as the kernels read them (valid while laid)
  [[nodiscard]] LayoutView view();

  /// the container grew from `old` slots (products placed, arrivals): the new slots are tail (-1)
  void grown(const pic::ParticleContainer& con, std::size_t old);

  /// widen the bounds of every key by every live LP it holds, on the caller's queue (the
  /// caller waits)
  void touch_all(pic::ParticleContainer& con, const tyvi::mdgrid_work& w);

  /// number of live LPs and their summed weight, energy and momentum
  struct Totals {
    double N, W, E;           // live LPs, their weight and energy
    std::array<double, 3> P;  // their momentum
  };
  [[nodiscard]] static Totals totals(const pic::ParticleContainer& con);

  // ---- host copies for IO and tests
  [[nodiscard]] std::vector<std::int32_t> get_slot_keys(const pic::ParticleContainer& con);
  struct Tables {
    std::vector<runko::index_t> b, cap, fill, prod;  // per key, as on the device
    std::vector<float> wmax, emin, emax;             // per key
    double llo, dl;                                  // log bin edges
  };
  [[nodiscard]] Tables get_tables();
};

}  // namespace qed
