#pragma once

// Due to the templated interpolator this can not be in its own compilation unit,
// i.e. this has to be in a header. To not make particle.h too long
// the pushers are in own separate headers which are included at the end of particle.h.

#include "core/pic/particle.h"
#include "tools/math.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

void
  pic::ParticleContainer::push_particles_higuera_cary(
    const double cfl_d,
    const runko::EB_interpolator<value_type> auto interpolator)
{
  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // particle (four-)velocities
  const auto ids_mds = ids_.mds();

  using vt = value_type;

  const vt cfl   = static_cast<vt>(cfl_d);          // c dx/dt
  const vt qm    = toolbox::sign(charge_) / mass_;  // charge-to-mass ratio
  const vt hqm   = vt { 0.5 } * qm;                 // half charge-to-mass ratio
  const vt cfl2  = cfl * cfl;                       // c^2
  const vt cinv  = vt { 1 } / cfl;                  // 1/c
  const vt cinv2 = cinv * cinv;                     // 1/c^2

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        if(ids_mds[idx][] == runko::dead_prtc_id) { return; }

        using Vec3 = toolbox::Vec3<vt>;

        const auto eb = interpolator(Vec3(pos_mds[idx]));
        const Vec3& E = eb.E;
        const Vec3& B = eb.B;

        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);
        const Vec3 E0 = hqm * E;
        const Vec3 u0 = v0 + E0;

        // B half-impulse (NOT divided by cfl — needed for ginv formula)
        const Vec3 Bt = hqm * B;

        // Higuera-Cary intermediate Lorentz factor
        const vt u0sq  = toolbox::dot(u0, u0);
        const vt b2    = toolbox::dot(Bt, Bt);
        const vt bdotu = toolbox::dot(Bt, u0);
        const vt gmb   = vt { 1 } + u0sq * cinv2 - b2 * cinv2;
        const vt disc  = gmb * gmb + vt { 4 } * (b2 * cinv2 + bdotu * bdotu * cinv2);
        const vt ginv  = vt { 1 } / sstd::sqrt(vt { 0.5 } * (gmb + sstd::sqrt(disc)));

        // Scale B by ginv/c to get rotation parameter
        const vt gc   = ginv * cinv;
        const Vec3 B0 = gc * Bt;

        const vt f = vt { 2 } / (vt { 1 } + gc * gc * b2);

        const Vec3 u1 = f * (u0 + toolbox::cross(u0, B0));
        const Vec3 u2 = u0 + toolbox::cross(u1, B0) + E0;

        const vt ginv2 = cfl / sstd::sqrt(cfl2 + toolbox::dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          vel_mds[idx][i] = u2[i] * cinv;
          pos_mds[idx][i] += u2[i] * ginv2;
        }
      })
    .wait();
}
