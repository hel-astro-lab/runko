// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Due to the templated interpolator this can not be in its own compilation unit,
// i.e. this has to be in a header. To not make particle.h too long
// the pushers are in own separate headers which are included at the end of particle.h.
//
// Radiative pushers: Higuera-Cary Lorentz push (copied from particle_higuera_cary.h),
// then a radiative stage at the mid-step velocity (Tamburini+10 splitting) with the
// gamma^2 term implicit in u, then an Euler-Maruyama kick for the stochastic closures.
// Formulas and notation follow radiative_drag_v2.tex.

#include "runko/pic/particle.h"
#include "runko/tools/math.h"
#include "runko/tools/rng.h"
#include "runko/tools/vector.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>

template<pic::ParticlePusher rad>
inline void
  pic::ParticleContainer::push_particles_radiation(
    const double cfl_d,
    const runko::EB_interpolator<value_type> auto interpolator,
    const RadParams params,
    const std::uint64_t cntr,
    const std::uint32_t species,
    const std::uint32_t seed)
{
  using enum ParticlePusher;
  constexpr bool is_sync = rad == rad_sync or rad == rad_sync_ssa;
  constexpr bool beam    = rad == rad_beam;
  constexpr bool heat    = rad == rad_comp_heat;
  constexpr bool ssa     = rad == rad_sync_ssa;

  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // particle (four-)velocities
  const auto ids_mds = ids_.mds();

  using vt   = value_type;
  using Vec3 = toolbox::Vec3<vt>;

  const vt cfl   = static_cast<vt>(cfl_d);                      // c dx/dt
  const vt qm    = static_cast<vt>(sstd::sign(charge_) / mass_);  // charge-to-mass ratio
  const vt hqm   = vt { 0.5 } * qm;                             // half charge-to-mass ratio
  const vt cfl2  = cfl * cfl;                                   // c^2
  const vt cinv  = vt { 1 } / cfl;                              // 1/c
  const vt cinv2 = cinv * cinv;                                 // 1/c^2

  // sigma_T/m ~ m^-3 so heavier species radiate less
  const double m3 = mass_ * mass_ * mass_;
  const vt A      = static_cast<vt>(params.drag / m3);
  // beam: 3/4 A_b dt and unit vector n from the vector A_b dt n
  const auto& b   = params.beam;
  const double bb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
  [[maybe_unused]] const vt Ab = static_cast<vt>(0.75 * bb / m3);
  [[maybe_unused]] const Vec3 nb =
    bb > 0.0 ? Vec3 { static_cast<vt>(b[0] / bb), static_cast<vt>(b[1] / bb), static_cast<vt>(b[2] / bb) }
             : Vec3 { vt { 0 }, vt { 0 }, vt { 0 } };
  // rad_temp is in m_e c^2; in units of this species' m c^2 it is rad_temp m_e/m
  const double theta_d = params.rad_temp / mass_;
  [[maybe_unused]] const vt theta = static_cast<vt>(theta_d);
  // SSA escape function tau(gamma) = (gamma_t/gamma)^{13/3} exp[-c_M (gamma^{2/3} - gamma_t^{2/3})]
  [[maybe_unused]] const vt lgt  = static_cast<vt>(std::log(params.gamma_t));
  [[maybe_unused]] const vt gt23 = static_cast<vt>(std::pow(params.gamma_t, 2.0 / 3.0));
  [[maybe_unused]] const vt cM   = static_cast<vt>(
    1.8899 * std::cbrt(std::numbers::pi / 4.0) * std::pow(std::max(theta_d, 1e-30), -2.0 / 3.0));
  [[maybe_unused]] constexpr vt eps = vt { 1e-12 };  // guards divisions by 0

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        if(ids_mds[idx][] == runko::dead_prtc_id) { return; }

        const auto eb = interpolator(Vec3(pos_mds[idx]));
        const Vec3& E = eb.E;
        const Vec3& B = eb.B;

        // --- Higuera-Cary Lorentz push (copy of particle_higuera_cary.h) ---
        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);
        const Vec3 E0 = hqm * E;
        const Vec3 u0 = v0 + E0;
        const Vec3 Bt = hqm * B;  // B half-impulse (NOT divided by cfl)

        const vt u0sq  = toolbox::dot(u0, u0);
        const vt b2    = toolbox::dot(Bt, Bt);
        const vt bdotu = toolbox::dot(Bt, u0);
        const vt gmb   = vt { 1 } + u0sq * cinv2 - b2 * cinv2;
        const vt disc  = gmb * gmb + vt { 4 } * (b2 * cinv2 + bdotu * bdotu * cinv2);
        const vt ginv  = vt { 1 } / sstd::sqrt(vt { 0.5 } * (gmb + sstd::sqrt(disc)));

        const vt gc   = ginv * cinv;
        const Vec3 B0 = gc * Bt;
        const vt f    = vt { 2 } / (vt { 1 } + gc * gc * b2);
        const Vec3 u1 = f * (u0 + toolbox::cross(u0, B0));
        const Vec3 uL = u0 + toolbox::cross(u1, B0) + E0;  // c u_L

        // --- radiative drift at the mid-step velocity u^n ---
        const Vec3 un = vt { 0.5 } * cinv * (uL + v0);
        const vt gn   = sstd::sqrt(vt { 1 } + toolbox::dot(un, un));
        const Vec3 bn = un / gn;  // beta^n

        vt kappa   = A * gn;  // compton: -A gamma^2 beta = -(A gamma) u
        Vec3 a_mid = Vec3 { vt { 0 }, vt { 0 }, vt { 0 } };  // c a_mid dt
        [[maybe_unused]] vt F2 = vt { 0 };  // (E + beta x B)^2 - (beta.E)^2; gamma^2 F2 = E_rest^2
        if constexpr(is_sync) {
          // Landau-Lifshitz reduced force (Vranic+16 eq. 9) in units of B_0
          const Vec3 fL = E + toolbox::cross(bn, B);
          const vt bE   = toolbox::dot(bn, E);
          F2            = toolbox::dot(fL, fL) - bE * bE;
          kappa *= F2;
          a_mid = (A * cfl) * (toolbox::cross(E, B)
                               + toolbox::cross(B, toolbox::cross(B, bn)) + bE * E);
        }
        if constexpr(beam) {
          // Compton rocket: pressure 3/4 A_b (1 - beta.n) n and drag 3/4 A_b gamma (1 - beta.n)^2
          const vt w = vt { 1 } - toolbox::dot(bn, nb);
          kappa += Ab * gn * w * w;
          a_mid = a_mid + (Ab * cfl * w) * nb;
        }

        // SSA: drift frame (E' || B'), boosted u', escape function, heating multiplier h
        [[maybe_unused]] vt betaD = 0, gammaD = 1, s = vt { 0 };
        [[maybe_unused]] Vec3 nD = Vec3 { vt { 1 }, vt { 0 }, vt { 0 } };
        if constexpr(ssa) {
          const Vec3 ExB = toolbox::cross(E, B);
          const vt S2    = toolbox::dot(ExB, ExB);
          const vt E2    = toolbox::dot(E, E);
          const vt B2    = toolbox::dot(B, B);
          const vt S     = sstd::sqrt(S2);
          const vt W     = E2 + B2;
          // v_D/c = [W - sqrt(W^2 - 4 S^2)] / 2S written without cancellation
          betaD  = vt { 2 } * S / (W + sstd::sqrt(sstd::max(W * W - vt { 4 } * S2, vt { 0 })) + eps);
          gammaD = vt { 1 } / sstd::sqrt(sstd::max(vt { 1 } - betaD * betaD, eps));  // finite as E -> B
          nD     = ExB / (S + eps);
          const vt upar = toolbox::dot(un, nD);
          const vt gp   = gammaD * (gn - betaD * upar);
          const vt pp2  = sstd::max(gp * gp - vt { 1 }, vt { 0 });
          // escape s = 1 - e^-tau and s' = ds/dgamma'; log tau clamped so tau e^-tau stays finite
          const vt lg   = sstd::log(gp);
          const vt g23  = sstd::exp(vt { 2 } / vt { 3 } * lg);
          const vt ltau = sstd::min(vt { 13 } / vt { 3 } * (lgt - lg) - cM * (g23 - gt23), vt { 30 });
          const vt tau  = sstd::exp(ltau);
          const vt etau = sstd::exp(-tau);
          const vt thick = E2 < B2 ? vt { 1 } : vt { 0 };  // no absorption where E > B
          s             = thick * (vt { 1 } - etau);
          const vt sp   = -thick * tau * etau
                        * (vt { 13 } / (vt { 3 } * gp) + vt { 2 } * cM / (vt { 3 } * sstd::exp(lg / vt { 3 })));
          // h of Eq. (h): 1 where emission escapes, < 0 where absorption heats; scales the whole LL force
          const vt h = vt { 1 } - theta * ((vt { 1 } / gp + vt { 3 } * gp / (pp2 + eps)) * s + sp);
          kappa *= h;
          a_mid = h * a_mid;
        }

        // implicit factor for cooling, explicit update for heating (kappa < 0)
        Vec3 u2 = kappa >= vt { 0 } ? (uL + a_mid) / (vt { 1 } + kappa)
                                    : (uL + a_mid) * (vt { 1 } - kappa);

        // --- stochastic kick, coefficients at u^n ---
        if constexpr(heat or ssa) {
          const auto xi = toolbox::normal4<vt>(ids_mds[idx][], cntr, species, seed);
          if constexpr(heat) {
            // Thomson diffusion tensor D_par bb + D_perp (1 - bb); gamma^2 beta^2 = p^2
            const vt AT    = A * theta;
            const vt gn2   = gn * gn;
            const vt p2    = toolbox::dot(un, un);
            const vt Dpar  = AT * gn2 * (gn2 + vt { 3.2 } * p2);
            const vt Dperp = AT * (gn2 - vt { 0.1 } * p2);
            // kick std capped at gamma^n: beyond 2 D dt ~ gamma^2 (8 A Theta gamma^4 > gamma^2)
            // the step is invalid and gamma^4 runs away to inf/NaN in a few laps
            const vt kpar  = sstd::min(sstd::sqrt(vt { 2 } * Dpar), gn);
            const vt kperp = sstd::min(sstd::sqrt(vt { 2 } * Dperp), gn);
            // D^{1/2} xi = kperp xi + (kpar - kperp)(bh.xi) bh; at rest kpar = kperp
            const Vec3 x  = Vec3 { xi[0], xi[1], xi[2] };
            const Vec3 bh = un / (sstd::sqrt(p2) + eps);
            u2 = u2 + cfl * (kperp * x + ((kpar - kperp) * toolbox::dot(bh, x)) * bh);
          }
          if constexpr(ssa) {
            // Eq. (kick_ssa): step gamma' in the drift frame, reflect at 1, rescale u' at fixed pitch;
            // D_e = Theta s A gamma^2 F2 (Eq. De_code). A du' = dgamma'/beta' step would add a spurious Ito drift
            const vt De    = theta * s * A * gn * gn * F2;
            const Vec3 u   = cinv * u2;
            const vt g     = sstd::sqrt(vt { 1 } + toolbox::dot(u, u));
            const vt upar  = toolbox::dot(u, nD);
            const vt gpk   = gammaD * (g - betaD * upar);
            const Vec3 up  = u + (gammaD * (upar - betaD * g) - upar) * nD;
            const vt ppk   = sstd::sqrt(toolbox::dot(up, up));
            const vt gpn   = vt { 1 } + sstd::abs(gpk - vt { 1 } + xi[0] * sstd::sqrt(vt { 2 } * De / gammaD));
            const vt ppn   = sstd::sqrt(gpn * gpn - vt { 1 });
            const Vec3 upn = (ppn / (ppk + eps)) * up;
            const vt uparn = toolbox::dot(upn, nD);
            u2 = cfl * (upn + (gammaD * (uparn + betaD * gpn) - uparn) * nD);
          }
        }

        const vt ginv2 = cfl / sstd::sqrt(cfl2 + toolbox::dot(u2, u2));
        for(auto i = 0uz; i < 3uz; ++i) {
          vel_mds[idx][i] = u2[i] * cinv;
          pos_mds[idx][i] += u2[i] * ginv2;
        }
      })
    .wait();
}
