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

#include <array>
#include <cstdint>
#include <vector>

namespace qed {

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
  float* emin;
  float* emax;
  int nbin;
  int ncell_edge;
  int ncell;
  int nkeys;
  bool keyed;
  float mass;
  double llo;           // log of the lowest kinetic energy of the epoch
  double dl;            // log width of an energy bin
  std::array<float, 3> mins;  // tile extent, for the reaction cell of a position
  std::array<float, 3> maxs;

  [[nodiscard]] constexpr int ncells() const { return ncell; }

  /// reaction cell of a position
  [[nodiscard]] constexpr int cell_of(const float x, const float y, const float z) const
  {
    const auto at = [&](const float v, const int i) {
      const auto f = (v - mins[i]) / (maxs[i] - mins[i]) * static_cast<float>(ncell_edge);
      return sstd::clamp(static_cast<int>(sstd::floor(f)), 0, ncell_edge - 1);
    };
    return (at(x, 0) * ncell_edge + at(y, 1)) * ncell_edge + at(z, 2);
  }

  /// energy bin of an LP of energy e: log-spaced in the kinetic energy of the epoch
  [[nodiscard]] constexpr int bin_of(const double e) const
  {
    if(nbin == 1) { return 0; }
    const auto v = sstd::max(e - static_cast<double>(mass), 1e-30);
    return sstd::clamp(static_cast<int>(sstd::floor((sstd::log(v) - llo) / dl)), 0, nbin - 1);
  }

  [[nodiscard]] constexpr int key_of(const int cell, const int ebin, const int bit) const
  {
    const auto k = cell * nbin + ebin;
    return keyed ? 2 * k + bit : k;
  }

  /// cell and bin of a key
  [[nodiscard]] constexpr int cell_of_key(const int k) const
  {
    return (keyed ? k / 2 : k) / nbin;
  }
  [[nodiscard]] constexpr int bin_of_key(const int k) const
  {
    return (keyed ? k / 2 : k) % nbin;
  }
};

class [[nodiscard]] Pop {
  // The layout of one species; the container itself belongs to the tile, which passes it to
  // every call, so that a tile may be copied or moved without a pointer going stale.
  bool keyed_;
  float mass_;
  LayoutParams par_ {};
  std::array<float, 3> mins_ {}, maxs_ {};
  double llo_ { 0.0 }, dl_ { 1.0 };
  int nkeys_ { 0 };
  bool laid_ { false };

  thrust::device_vector<std::int32_t> skey_;
  thrust::device_vector<runko::index_t> b_, cap_, fill_, prod_;
  thrust::device_vector<float> wmax_, emin_, emax_;
  // scratch
  thrust::device_vector<runko::index_t> keys_, ord_, cnt_, slots_;

  /// widen the bounds of every key by every live LP it holds
  void touch_all(pic::ParticleContainer& con);

public:
  Pop(bool keyed, double mass);

  [[nodiscard]] bool laid() const { return laid_; }
  [[nodiscard]] int nkeys() const { return nkeys_; }
  [[nodiscard]] bool keyed() const { return keyed_; }
  [[nodiscard]] double mass() const { return mass_; }
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

  /// keep the key array in step with a container that grew: new slots are tail (-1)
  void follow_container(const pic::ParticleContainer& con);

  /// number of live LPs and their summed weight, energy and momentum
  struct Totals {
    double N, W, E;
    std::array<double, 3> P;
  };
  [[nodiscard]] static Totals totals(const pic::ParticleContainer& con);

  // ---- host copies for IO and tests
  [[nodiscard]] std::vector<std::int32_t> get_slot_keys(const pic::ParticleContainer& con);
  struct Tables {
    std::vector<runko::index_t> b, cap, fill, prod;
    std::vector<float> wmax, emin, emax;
    double llo, dl;
  };
  [[nodiscard]] Tables get_tables();
};

}  // namespace qed
