// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// Radiative pushers: Higuera-Cary Lorentz push followed by the cooling drift at the 
// mid-step velocity (Tamburini+10 splitting). For rad_comp_heat / rad_sync_ssa a 
// fluctuation-dissipation kick whose diffusion follows from the drift by the Einstein relation. 

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

  using toolbox::cross;
  using toolbox::dot;

  constexpr bool is_sync = rad == rad_sync or rad == rad_sync_ssa;
  constexpr bool beam    = rad == rad_beam;
  constexpr bool heat    = rad == rad_comp_heat;
  constexpr bool ssa     = rad == rad_sync_ssa;

  const auto pos_mds = pos_.mds();  // particle positions
  const auto vel_mds = vel_.mds();  // particle (four-)velocities
  const auto ids_mds = ids_.mds();

  using vt   = value_type;
  using Vec3 = toolbox::Vec3<vt>;

  const vt cfl   = static_cast<vt>(cfl_d);                        // c dx/dt
  const vt qm    = static_cast<vt>(sstd::sign(charge_) / mass_);  // charge-to-mass ratio
  const vt hqm   = vt { 0.5 } * qm;                               // half charge-to-mass ratio
  const vt cfl2  = cfl * cfl;                                     // c^2
  const vt cinv  = vt { 1 } / cfl;                                // 1/c
  const vt cinv2 = cinv * cinv;                                   // 1/c^2

  // sigma_T/m ~ m^-3 so heavier species radiate less
  const double m3 = mass_ * mass_ * mass_;
  const vt A      = static_cast<vt>(params.drag / m3);
  const auto& b   = params.beam; // beam: 3/4 A_b dt and unit vector n from the vector A_b dt n
  const double bb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
  [[maybe_unused]] const vt Ab = static_cast<vt>(0.75 * bb / m3);
  [[maybe_unused]] const Vec3 nb =
    bb > 0.0 ? Vec3 { static_cast<vt>(b[0] / bb), static_cast<vt>(b[1] / bb), static_cast<vt>(b[2] / bb) }
             : Vec3 { vt { 0 }, vt { 0 }, vt { 0 } };
  const double theta_d = params.rad_temp / mass_; // rad_temp in m_e c^2
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

        // --- Higuera-Cary Lorentz push ---
        const Vec3 v0 = cfl * Vec3(vel_mds[idx]);
        const Vec3 E0 = hqm * E;
        const Vec3 u0 = v0 + E0;
        const Vec3 Bt = hqm * B;  // B half-impulse (NOT divided by cfl)

        const vt u0sq  = dot(u0, u0);
        const vt b2    = dot(Bt, Bt);
        const vt bdotu = dot(Bt, u0);
        const vt gmb   = vt { 1 } + u0sq * cinv2 - b2 * cinv2;
        const vt disc  = gmb * gmb + vt { 4 } * (b2 * cinv2 + bdotu * bdotu * cinv2);
        const vt ginv  = vt { 1 } / sstd::sqrt(vt { 0.5 } * (gmb + sstd::sqrt(disc)));

        const vt gc   = ginv * cinv;
        const Vec3 B0 = gc * Bt;
        const vt f    = vt { 2 } / (vt { 1 } + gc * gc * b2);
        const Vec3 u1 = f * (u0 + cross(u0, B0));
        const Vec3 uL = u0 + cross(u1, B0) + E0;  // c u_L

        // --- radiative drift at the mid-step velocity u^n ---
        const Vec3 un = vt { 0.5 } * cinv * (uL + v0);
        const vt gn   = sstd::sqrt(vt { 1 } + dot(un, un));
        const Vec3 bn = un / gn;  // beta^n

        vt kappa   = A * gn;  // compton: -A gamma^2 beta = -(A gamma) u
        Vec3 a_mid = Vec3 { vt { 0 }, vt { 0 }, vt { 0 } };  // c a_mid dt

        //----------------------------------------------------------------------------------------- 
        [[maybe_unused]] vt F2 = vt { 0 };  // (E + beta x B)^2 - (beta.E)^2; gamma^2 F2 = E_rest^2
        [[maybe_unused]] vt kB = vt { 0 };  // B x (B x beta) = -(B^2/gamma) u_perp: friction across B
        [[maybe_unused]] Vec3 bh = Vec3 { vt { 0 }, vt { 0 }, vt { 0 } };  // B/|B|
        if constexpr(is_sync) { // Landau-Lifshitz reduced force (Vranic+16 eq. 9) in units of B_0
          const Vec3 fL = E + cross(bn, B);
          const vt bE   = dot(bn, E);
          const vt B2   = dot(B, B);
          F2            = dot(fL, fL) - bE * bE;
          kappa *= F2;
          kB    = A * B2 / gn;
          bh    = B / (sstd::sqrt(B2) + eps);
          a_mid = (A * cfl) * (cross(E, B) + bE * E);
        }

        //----------------------------------------------------------------------------------------- 
        if constexpr(beam) { // Beamed radiation field; Compton rocket
          const vt w = vt { 1 } - dot(bn, nb);
          kappa += Ab * gn * w * w;
          a_mid = a_mid + (Ab * cfl * w) * nb;
        }

        //----------------------------------------------------------------------------------------- 
        // cooling: every friction implicit, u2 = (1 + K)^{-1} (uL + a_mid), K = kappa + kB (1 - bh bh)
        const Vec3 uw = uL + a_mid;
        Vec3 u2       = uw / (vt { 1 } + kappa + kB);
        if constexpr(is_sync) {
          u2 = u2 + (dot(uw, bh) * (vt { 1 } / (vt { 1 } + kappa) - vt { 1 } / (vt { 1 } + kappa + kB))) * bh;
        }

        //----------------------------------------------------------------------------------------- 
        // Heating: Euler-Maruyama kick, coefficients at u^n ---
        if constexpr(heat or ssa) {
          const auto xi = toolbox::normal4<vt>(ids_mds[idx][], cntr, species, seed);
          const Vec3 x  = Vec3 { xi[0], xi[1], xi[2] };

          //----------------------------------------------------------------------------------------- 
          if constexpr(heat) { // Compton heating; isotropic, box frame
            const vt D = A * theta * (gn * gn + vt { 2 } * theta * gn + vt { 2 } * theta * theta);
            u2 = u2 + (cfl * sstd::sqrt(vt { 2 } * D)) * x;
          }

          //----------------------------------------------------------------------------------------- 
          if constexpr(ssa) { // SSA: drift frame (E' || B'); boost along nD
            const Vec3 ExB = cross(E, B);
            const vt E2    = dot(E, E);
            const vt B2    = dot(B, B);
            const vt EB    = dot(E, B);
            const vt S     = sstd::sqrt(dot(ExB, ExB));
            const vt W     = E2 + B2;
            // v_D/c = [W - sqrt(W^2 - 4 S^2)] / 2S written without cancellation
            const vt betaD  = vt { 2 } * S / (W + sstd::sqrt(sstd::max(W * W - vt { 4 } * S * S, vt { 0 })) + eps);
            const vt gammaD = vt { 1 } / sstd::sqrt(sstd::max(vt { 1 } - betaD * betaD, eps));
            const Vec3 nD   = ExB / (S + eps);
            const vt Bp2    = vt { 0.5 } * ((B2 - E2) + sstd::sqrt((B2 - E2) * (B2 - E2) + vt { 4 } * EB * EB));

            // lorentz boost helper 
            const auto boost = [=](const Vec3& u, const vt sgn) {  // sgn = -1: to drift frame, +1: back
              const vt g  = sstd::sqrt(vt { 1 } + dot(u, u));
              const vt up = dot(u, nD);
              return u + (gammaD * (up + sgn * betaD * g) - up) * nD;
            };

            // escape s(gamma'); log tau clamped; off where E' > B' (= E > B)
            const vt gp    = gammaD * (gn - betaD * dot(un, nD));
            const vt lg    = sstd::log(gp);
            const vt ltau  = sstd::min(vt { 13 } / vt { 3 } * (lgt - lg) - cM * (sstd::exp(vt { 2 } / vt { 3 } * lg) - gt23), vt { 30 });
            const vt thick = E2 < B2 ? vt { 1 } : vt { 0 };
            const vt s     = thick * (vt { 1 } - sstd::exp(-sstd::exp(ltau)));

            // D_s,perp and D_s,par with the invariant gamma^2 F2, times s, per box time (1/gammaD)
            const vt g2F2  = gn * gn * F2;
            const vt w     = vt { 2 } * A * theta * s / gammaD;
            const vt kpar  = sstd::sqrt(w * g2F2);
            const vt kperp = sstd::sqrt(w * (g2F2 + Bp2 * (vt { 1 } + vt { 2 } * theta * gp + vt { 2 } * theta * theta)));

            // D^{1/2} x = k_perp x + (k_par - k_perp)(x.bh) bh
            const Vec3 du = kperp * x + ((kpar - kperp) * dot(x, bh)) * bh;
            u2 = cfl * boost(boost(cinv * u2, vt { -1 }) + du, vt { 1 });
          }
        }

        //----------------------------------------------------------------------------------------- 
        const vt ginv2 = cfl / sstd::sqrt(cfl2 + dot(u2, u2));
        for(auto i = 0uz; i < 3uz; ++i) {
          vel_mds[idx][i] = u2[i] * cinv;
          pos_mds[idx][i] += u2[i] * ginv2;
        }
      })
    .wait();
}
