#include "core/pic/particle.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

#include <cmath>

void
  pic::ParticleContainer::push_particles_boris(
    const value_type cfl,
    const InterpolatedEB_function& EBfunc)
{
  const auto EB = EBfunc(pos_);

  const auto pos_mds = pos_.mds();
  const auto vel_mds = vel_.mds();
  const auto Emds    = EB.E.mds();
  const auto Bmds    = EB.B.mds();

  const auto qm = toolbox::sign(charge_) / mass_;

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        using Vec3 = toolbox::Vec3<value_type>;

        const value_type v0 = cfl * Vec3(vel_mds[idx]);
        const value_type E0 = 0.5f * qm * Vec3(Emds[idx]);
        const value_type u0 = v0 + E0;

        const value_type ginv = cfl / std::sqrt(cfl * cfl + toolbox::dot(u0, u0));

        const value_type B0 = 0.5f * qm * ginv * Vec3(Bmds[idx]) / cfl;

        const value_type f = 2.f / (1.f + toolbox::dot(B0, B0));

        const value_type u1 = f * (u0 + toolbox::cross(u0, B0));
        const value_type u2 = u0 + toolbox::cross(u1, B0) + E0;

        for(auto i = 0uz; i < 3uz; ++i) { vel_mds[idx][i] = u2[i] / cfl; }

        const value_type ginv2 = cfl / std::sqrt(cfl * cfl + toolbox::dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
