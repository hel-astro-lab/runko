#include "core/pic2/particle.h"
#include "tools/signum.h"
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
        using vec3 = std::array<double, 3>;

        auto v0 = vec3 {};
        auto E0 = vec3 {};
        auto u0 = vec3 {};
        auto B0 = vec3 {};
        auto u1 = vec3 {};
        auto u2 = vec3 {};

        const auto dot2 = [](const auto& m) {
          return m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
        };

        const auto cross = [](const auto& a, const auto& b) {
          return std::array<double, 3> { a[1] * b[2] - a[2] * b[1],
                                         a[0] * b[2] - a[2] * b[0],
                                         a[0] * b[1] - a[1] * b[0] };
        };

        for(auto i = 0uz; i < 3uz; ++i) {
          v0[i] = cfl * vel_mds[idx][i];
          E0[i] = 0.5 * qm * Emds[idx][i];
          u0[i] = v0[i] + E0[i];
        }

        const auto ginv = cfl / std::sqrt(cfl * cfl + dot2(u0));

        for(auto i = 0uz; i < 3uz; ++i) {
          B0[i] = 0.5 * qm * Bmds[idx][i] * ginv / cfl;
        }

        const auto f = 2 / (1 + dot2(B0));

        const auto u0crossB0 = cross(u0, B0);
        for(auto i = 0uz; i < 3uz; ++i) { u1[i] = f * (u0[i] + u0crossB0[i]); }

        const auto u1crossB0 = cross(u1, B0);
        for(auto i = 0uz; i < 3uz; ++i) {
          u2[i]           = u0[i] + u1crossB0[i] + E0[i];
          vel_mds[idx][i] = u2[i] / cfl;
        }

        const auto ginv2 = cfl / std::sqrt(cfl * cfl + dot2(u0));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
