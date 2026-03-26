// Faraday-Cayley particle pusher
//
// Applies the full electromagnetic bivector F = E + IB as a single Lorentz
// transformation via the Cayley rational parametrization, eliminating the
// E/B Strang splitting error of Boris/Higuera-Cary.
//
// Algorithm (Sherman-Morrison reduction of the 4x4 Cayley to 3x3):
//   1. Half E-kick to estimate gamma_eff (same as Boris)
//   2. Compute dimensionless parameters eps = kappa*E, beta = kappa*B
//      where kappa = (q dt)/(2 m gamma_eff c) = 0.5*qm/geff_cfl
//   3. Build RHS: w0 = gcfl + eps.v0, w = v0 + eps*gcfl + v0 x beta
//   4. Boris-like inverse of (I + [beta x]) applied to W = w + w0*eps
//   5. Sherman-Morrison correction for the rank-1 boost perturbation
//   6. Position update (identical to Boris)
//
// Properties:
//   - Exact Lorentz transformation (mass shell preserved to machine precision)
//   - Second-order accurate (same order as Boris/HC)
//   - Reduces to Boris exactly when E = 0
//   - No E/B splitting error (Thomas precession handled correctly)
//   - Branchless (no case analysis for field geometry)

#include "core/pic/particle.h"
#include "tools/signum.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

#include "tools/math.h"

void
  pic::ParticleContainer::push_particles_faraday(
    const double cfl_,
    const InterpolatedEB_function& EBfunc)
{
  using vt = value_type;
  using toolbox::dot;
  using toolbox::cross;

  const auto EB = EBfunc(pos_);

  const auto pos_mds = pos_.mds();
  const auto vel_mds = vel_.mds();
  const auto Emds    = EB.E.mds();
  const auto Bmds    = EB.B.mds();

  const vt cfl = cfl_;
  const vt qm  = toolbox::sign(charge_) / mass_;

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        using Vec3 = toolbox::Vec3<vt>;

        const Vec3 E = Vec3(Emds[idx]);
        const Vec3 B = Vec3(Bmds[idx]);

        // current four-velocity in code units
        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);

        // gamma * cfl
        const vt gcfl = sstd::sqrt(cfl * cfl + dot(v0, v0));

        // half E-kick for gamma_eff (Boris prescription)
        const Vec3 u0       = v0 + vt { 0.5 } * qm * E;
        const vt   geff_cfl = sstd::sqrt(cfl * cfl + dot(u0, u0));

        // Cayley parameters: kappa = (q dt) / (2 m gamma_eff c)
        const vt   kappa = vt { 0.5 } * qm / geff_cfl;
        const Vec3 eps   = kappa * E;
        const Vec3 beta  = kappa * B;

        // RHS: (I + Omega/2) applied to v4 = (gcfl, v0)
        const vt   w0 = gcfl + dot(eps, v0);
        const Vec3 W  = v0 + eps * gcfl + cross(v0, beta) + w0 * eps;

        // Boris-like inverse: (I + [beta x])^{-1} x
        //   = (x - beta x x + (beta.x) beta) / (1 + |beta|^2)
        const vt   b2 = dot(beta, beta);
        const vt   f  = vt { 1 } / (vt { 1 } + b2);

        const Vec3 W_rot   = f * (W - cross(beta, W) + dot(beta, W)*beta);

        // Sherman-Morrison correction for the rank-1 boost perturbation
        const vt   bde     = dot(beta, eps);
        const Vec3 eps_rot = f * (eps - cross(beta, eps) + bde * beta);

        const vt   D  = vt { 1 } - f * (dot(eps, eps) + bde * bde);
        const Vec3 u2 = W_rot + eps_rot * (dot(eps, W_rot) / D);

        // store final velocity (u/c)
        for(auto i = 0uz; i < 3uz; ++i) { vel_mds[idx][i] = u2[i] / cfl; }

        // position update: x += (u / gamma) * dt
        const vt ginv2 = cfl / sstd::sqrt(cfl * cfl + dot(u2, u2));

        for(auto i = 0uz; i < 3uz; ++i) {
          pos_mds[idx][i] = pos_mds[idx][i] + vel_mds[idx][i] * ginv2 * cfl;
        }
      })
    .wait();
}
