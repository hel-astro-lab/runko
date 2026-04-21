#pragma once

// Due to the templated interpolator this can not be in its own compilation unit,
// i.e. this has to be in a header. To not make particle.h too long
// the pushers are in own separate headers which are included at the end of particle.h.

#include "core/pic/particle.h"
#include "tools/math.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

inline void
  pic::ParticleContainer::push_particles_boris(
    const double cfl_d,
    const runko::EB_interpolator<value_type> auto interpolator)
{
  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // particle (four-)velocities
  const auto ids_mds = ids_.mds();

  const auto cfl =
    static_cast<value_type>(cfl_d);  // c dx/dt = speed of light in num. units
  const value_type qm = toolbox::sign(charge_) / mass_;  // charge-to-mass ratio

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        if(ids_mds[idx][] == runko::dead_prtc_id) { return; }

        using Vec3 = toolbox::Vec3<value_type>;

        const auto EB = interpolator(Vec3(pos_mds[idx]));

        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);
        const Vec3 E0 = value_type { 0.5 } * qm * EB.E;
        const Vec3 u0 = v0 + E0;

        const value_type ginv = cfl / sstd::sqrt(cfl * cfl + toolbox::dot(u0, u0));

        const Vec3 B0 = value_type { 0.5 } * qm * ginv * EB.B / cfl;

        const value_type f =
          value_type { 2 } / (value_type { 1 } + toolbox::dot(B0, B0));

        const Vec3 u1 = f * (u0 + toolbox::cross(u0, B0));
        const Vec3 u2 = u0 + toolbox::cross(u1, B0) + E0;

        for(auto i = 0uz; i < 3uz; ++i) { vel_mds[idx][i] = u2[i] / cfl; }

        const value_type ginv2 = cfl / sstd::sqrt(cfl * cfl + toolbox::dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
