#include "core/pic/particle.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

#include <cmath>

void
  pic::ParticleContainer::push_particles_boris(
    const double cfl_,
    const InterpolatedEB_function& EBfunc)
{
  const auto EB = EBfunc(pos_);

  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // particle (four-)velocities
  const auto Emds    = EB.E.mds();  // electric field vector
  const auto Bmds    = EB.B.mds();  // magnetic field vector

  const value_type cfl =
    cfl_;  // c dx/dt = speed of light in num. units (note: implicit type cast)
  const value_type qm = toolbox::sign(charge_) / mass_;  // charge-to-mass ratio

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        using Vec3 = toolbox::Vec3<value_type>;

        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);
        const Vec3 E0 = value_type { 0.5 } * qm * Vec3(Emds[idx]);
        const Vec3 u0 = v0 + E0;

        const value_type ginv = cfl / std::sqrt(cfl * cfl + toolbox::dot(u0, u0));

        const Vec3 B0 = value_type { 0.5 } * qm * ginv * Vec3(Bmds[idx]) / cfl;

        const value_type f =
          value_type { 2 } / (value_type { 1 } + toolbox::dot(B0, B0));

        const Vec3 u1 = f * (u0 + toolbox::cross(u0, B0));
        const Vec3 u2 = u0 + toolbox::cross(u1, B0) + E0;

        for(auto i = 0uz; i < 3uz; ++i) { vel_mds[idx][i] = u2[i] / cfl; }

        const value_type ginv2 = cfl / std::sqrt(cfl * cfl + toolbox::dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
