// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The species of a tile as the sampler kernels see them: device views of the three
// containers (photons, electrons, positrons) with their layouts, and the tile around
// them.

#include "runko/pic/particle.h"
#include "runko/qed/field.h"
#include "runko/qed/layout.h"
#include "runko/tools/math.h"
#include "thrust/device_vector.h"

#include <array>
#include <cstddef>
#include <utility>

namespace qed {

using I = runko::index_t;

/// Device pointer to the storage of v, valid until v is resized.
template<typename T>
[[nodiscard]] T*
  ptr(thrust::device_vector<T>& v)
{ return thrust::raw_pointer_cast(v.data()); }

/// Species roles; each container is found from its charge and mass.
enum class Species : std::size_t { ph = 0, em = 1, ep = 2 };

using VecMds = decltype(std::declval<pic::ParticleContainer&>().vel_mds());
using IdMds  = decltype(std::declval<pic::ParticleContainer&>().ids_mds());
using WgtMds = decltype(weight_mds(std::declval<pic::ParticleContainer&>()));

/// One species as the kernels see it, captured by value.
struct PView {
  VecMds pos, vel;  // positions (cells) and momenta (m_e c; k for photons)
  IdMds ids;        // dead_prtc_id marks a free slot
  WgtMds wgt;       // weights (payload component 0), empty for leptons
  bool weighted;    // wgt (photons) or w0 (leptons)
  float w0;         // common weight of a species without payload
  float mass;       // 0 for photons, 1 for leptons
  LayoutView lay;   // key tables, for touch()

  [[nodiscard]] constexpr bool alive(const I i) const
  { return ids[i][] != runko::dead_prtc_id; }

  [[nodiscard]] constexpr Vec3 u(const I i) const
  {
    return Vec3(
      static_cast<double>(vel[i][0]),
      static_cast<double>(vel[i][1]),
      static_cast<double>(vel[i][2]));
  }

  [[nodiscard]] constexpr double w(const I i) const
  { return static_cast<double>(weighted ? wgt[i][] : w0); }

  [[nodiscard]] constexpr double e(const I i) const
  { return static_cast<double>(energy_of(mass, vel[i][0], vel[i][1], vel[i][2])); }

  constexpr void set_u(const I i, const Vec3& v) const
  {
    for(auto k = 0uz; k < 3uz; ++k) { vel[i][k] = static_cast<float>(v[k]); }
  }

  constexpr void set_w(const I i, const double v) const
  {
    if(weighted) { wgt[i][] = static_cast<float>(v); }
  }

  constexpr void kill(const I i) const
  {
    ids[i][] = runko::dead_prtc_id;
    if(weighted) { wgt[i][] = 0.0f; }
  }

  /// Widen the key bounds over a slot updated in place. A slot kicked out of its bin
  /// leaves the key, inert until the next layout, so a bin's bounds stay its own.
  constexpr void touch(const I i) const
  {
    const auto k = lay.skey[i];
    if(k < 0) { return; }

    const auto ee = static_cast<float>(e(i));
    if(lay.bin_of(ee) != lay.bin_of_key(k)) {
      lay.skey[i] = -1;
      return;
    }

    lay.widen(k, ee, static_cast<float>(w(i)));
  }
};

/// View of a container alone, for the sources, the sinks and the pair sum.
inline PView
  make_view(pic::ParticleContainer& con, const double w0)
{
  return PView { .pos      = con.pos_mds(),
                 .vel      = con.vel_mds(),
                 .ids      = con.ids_mds(),
                 .wgt      = weight_mds(con),
                 .weighted = con.n_payload() > 0uz,
                 .w0       = static_cast<float>(w0),
                 .mass     = static_cast<float>(con.args().mass),
                 .lay      = {} };
}

/// View of a container and its layout; stale after a resize.
inline PView
  make_view(pic::ParticleContainer& con, Layout& layout)
{
  auto v = make_view(con, layout.w0());
  v.lay  = layout.view();
  return v;
}

/// What the sampler needs of the tile for one step.
struct [[nodiscard]] TileView {
  FieldView field;
  std::array<pic::ParticleContainer*, 3> con;  // ph, e-, e+
  std::array<Layout*, 3> layout;                    // their layouts
  std::array<float, 3> mins, maxs;             // tile extent
  std::array<runko::prtc_id_type, 3> id_base;  // first free id per species
  std::array<std::size_t, 3> ids_used;         // returned: ids consumed by the step
};

/// Containers, layouts and device views of one step.
struct [[nodiscard]] Ctx {
  TileView* t;             // containers, layouts, extent and ids
  std::array<PView, 3> v;  // views of ph, e-, e+; stale after a resize
  int ncell;               // reaction cells of the tile
  int nb;                  // energy bins per cell

  Ctx(TileView& tile, const int nbin) : t { &tile }, v {}, ncell { 0 }, nb { nbin }
  {
    refresh_all();
    ncell = v[0].lay.ncell;
  }

  void refresh_all()
  {
    for(auto i = 0uz; i < 3uz; ++i) { v[i] = make_view(*t->con[i], *t->layout[i]); }
  }
};

}  // namespace qed
