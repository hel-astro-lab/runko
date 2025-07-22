#include "core/pic2/particle.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

#include <cmath>

void
  pic2::ParticleContainer::push_particles_boris(
    const double cfl,
    const InterpolatedEB_function EBfunc)
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
        using Vec3 = toolbox::Vec3<double>;

        const auto v0 = cfl * Vec3(vel_mds[idx]);
        const auto E0 = 0.5 * qm * Vec3(Emds[idx]);
        const auto u0 = v0 + E0;

        const auto ginv = cfl / std::sqrt(cfl * cfl + toolbox::dot(u0, u0));

        const auto B0 = 0.5 * qm * ginv * Vec3(Bmds[idx]) / cfl;

        const auto f = 2 / (1 + toolbox::dot(B0, B0));

        const auto u1 = f * (u0 + toolbox::cross(u0, B0));
        const auto u2 = u0 + toolbox::cross(u1, B0) + E0;

        for(auto i = 0uz; i < 3uz; ++i) { vel_mds[idx][i] = u2[i] / cfl; }

        const auto ginv2 = cfl / std::sqrt(cfl * cfl + toolbox::dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
