// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Included at the end of particle.h next to the other pushers.

#include "runko/pic/particle.h"
#include "runko/tools/math.h"
#include "runko/tools/vector.h"
#include "tyvi/mdgrid.h"

inline void
  pic::ParticleContainer::push_particles_photon(const double cfl_d)
{
  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // photon momenta k (energy x = |k|)
  const auto ids_mds = ids_.mds();

  const auto cfl = static_cast<value_type>(cfl_d);  // c dt/dx, the distance of a step

  // A photon moves at the speed of light along its momentum, whatever its energy:
  // pos += k/|k| c dt. Its momentum is untouched (no charge, no Lorentz force).
  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        if(ids_mds[idx][] == runko::dead_prtc_id) { return; }

        using Vec3    = toolbox::Vec3<value_type>;
        const Vec3 k  = Vec3(vel_mds[idx]);       // photon momentum
        const auto kk = toolbox::dot(k, k);       // |k|^2
        if(kk <= value_type { 0 }) { return; }    // no direction to move along
        const auto step = cfl / sstd::sqrt(kk);   // c dt per unit of k

        for(auto i = 0uz; i < 3uz; ++i) { pos_mds[idx][i] = pos_mds[idx][i] + k[i] * step; }
      })
    .wait();
}
