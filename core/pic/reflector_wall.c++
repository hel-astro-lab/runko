#include "core/pic/reflector_wall.h"

#include "core/emf/yee_lattice.h"
#include "core/pic/particle.h"
#include "core/pic/tile.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include "tyvi/sstd.h"

#include <algorithm>
#include <cstdint>
#include <ranges>

namespace {

using value_type = pic::ParticleContainer::value_type;
using Vec3       = toolbox::Vec3<value_type>;
using Vec3uz     = toolbox::Vec3<uint32_t>;
using toolbox::dot;

constexpr value_type EPS = 1e-10f;

}  // namespace


namespace pic {

/// Atomic zigzag current deposit for a single-particle sub-trajectory x1 → x2.
///
/// x1 and x2 are in lattice-local coordinates (relative to lattice origin).
/// Mirrors the atomic variant in particle_current_zigzag_1st.c++:276-329.
template<typename Jmds_t>
__attribute__((always_inline)) constexpr void
  zigzag_deposit_single(
    const Jmds_t Jmds,
    const Vec3 x1,
    const Vec3 x2,
    const value_type charge)
{
  const auto fi1 = Vec3(sstd::floor(x1(0)), sstd::floor(x1(1)), sstd::floor(x1(2)));
  const auto fi2 = Vec3(sstd::floor(x2(0)), sstd::floor(x2(1)), sstd::floor(x2(2)));

  // Relay point: branchless min/max (no lambda return → SIMD-safe)
  const auto relay = [&](const uint32_t j) {
    const auto a = sstd::min(fi1(j), fi2(j)) + value_type { 1 };
    const auto b =
      sstd::max(sstd::max(fi1(j), fi2(j)), value_type { 0.5 } * (x1(j) + x2(j)));
    return sstd::min(a, b);
  };

  const auto x_relay = Vec3(relay(0), relay(1), relay(2));

  const auto F1 = charge * (x_relay - x1);
  const auto F2 = charge * (x2 - x_relay);

  const auto i1 = fi1.template as<uint32_t>();
  const auto i2 = fi2.template as<uint32_t>();

  const auto W1 = value_type { 0.5 } * (x1 + x_relay) - fi1;
  const auto W2 = value_type { 0.5 } * (x2 + x_relay) - fi2;

  const auto [Fx1, Fy1, Fz1] = F1.data;
  const auto [Fx2, Fy2, Fz2] = F2.data;
  const auto [Wx1, Wy1, Wz1] = W1.data;
  const auto [Wx2, Wy2, Wz2] = W2.data;

  const auto store_current = [&](const Vec3uz& index, const Vec3& current) {
    const auto si  = index.template as<uint32_t>();
    auto* const Jx = &thrust::raw_reference_cast(Jmds[si.data][0]);
    auto* const Jy = &thrust::raw_reference_cast(Jmds[si.data][1]);
    auto* const Jz = &thrust::raw_reference_cast(Jmds[si.data][2]);

    sstd::atomic_add(Jx, current[0]);
    sstd::atomic_add(Jy, current[1]);
    sstd::atomic_add(Jz, current[2]);
  };

  static constexpr auto one = value_type { 1 };

  // clang-format off
  store_current(i1,                   Vec3(Fx1*(one - Wy1)*(one - Wz1),
                                           Fy1*(one - Wx1)*(one - Wz1),
                                           Fz1*(one - Wx1)*(one - Wy1)));
  store_current(i2,                   Vec3(Fx2*(one - Wy2)*(one - Wz2),
                                           Fy2*(one - Wx2)*(one - Wz2),
                                           Fz2*(one - Wx2)*(one - Wy2)));
  store_current(i1 + Vec3uz(1, 0, 0), Vec3(0,
                                           Fy1*Wx1*(one - Wz1),
                                           Fz1*Wx1*(one - Wy1)));
  store_current(i2 + Vec3uz(1, 0, 0), Vec3(0,
                                           Fy2*Wx2*(one - Wz2),
                                           Fz2*Wx2*(one - Wy2)));
  store_current(i1 + Vec3uz(0, 1, 0), Vec3(Fx1*Wy1*(one - Wz1),
                                           0,
                                           Fz1*(one - Wx1)*Wy1));
  store_current(i2 + Vec3uz(0, 1, 0), Vec3(Fx2*Wy2*(one - Wz2),
                                           0,
                                           Fz2*(one - Wx2)*Wy2));
  store_current(i1 + Vec3uz(0, 0, 1), Vec3(Fx1*(one - Wy1)*Wz1,
                                           Fy1*(one - Wx1)*Wz1,
                                           0));
  store_current(i2 + Vec3uz(0, 0, 1), Vec3(Fx2*(one - Wy2)*Wz2,
                                           Fy2*(one - Wx2)*Wz2,
                                           0));
  store_current(i1 + Vec3uz(0, 1, 1), Vec3(Fx1*Wy1*Wz1, 0, 0));
  store_current(i2 + Vec3uz(0, 1, 1), Vec3(Fx2*Wy2*Wz2, 0, 0));
  store_current(i1 + Vec3uz(1, 0, 1), Vec3(0, Fy1*Wx1*Wz1, 0));
  store_current(i2 + Vec3uz(1, 0, 1), Vec3(0, Fy2*Wx2*Wz2, 0));
  store_current(i1 + Vec3uz(1, 1, 0), Vec3(0, 0, Fz1*Wx1*Wy1));
  store_current(i2 + Vec3uz(1, 1, 0), Vec3(0, 0, Fz2*Wx2*Wy2));
  // clang-format on
}


/// Reflect particles behind the wall and deposit correction currents.
///
/// For each particle behind the wall: unwind to pre-push position, find the
/// collision point, reflect ux (relativistic Lorentz boost for moving walls),
/// and advance to the reflected position. Dual correction currents (+q real
/// path, -q cancellation) are deposited into correction_J so that no net
/// current remains behind the wall after the normal deposit pass.
void
  ParticleContainer::reflect_at_wall(
    const reflector_wall& wall,
    runko::VecGrid<emf::YeeLattice::value_type>& correction_J,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl)
{
  const auto pos_mds = pos_.mds();
  const auto vel_mds = vel_.mds();
  const auto Jmds    = correction_J.mds();

  const value_type walloc    = wall.walloc;
  const value_type betawall  = wall.betawall;
  const value_type gammawall = wall.gammawall;
  const value_type charge    = charge_;
  const value_type c         = cfl;

  // previous wall location
  const value_type walloc0 = walloc - betawall * c;

  const auto origo = Vec3(lattice_origo_coordinates);

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        // Branchless reflector wall: all particles execute the full
        // computation; float masks zero out effects for non-affected particles.
        //
        // Position convention:
        //   pos_1    = current post-push position (may be behind wall)
        //   pos_0    = unwound pre-push position (start of timestep)
        //   pos_col  = wall collision point
        //   pos_refl = final reflected position (in front of wall)
        const Vec3 pos_1        = Vec3(pos_mds[idx]);
        const Vec3 u            = Vec3(vel_mds[idx]);
        const value_type gam    = sstd::sqrt(1.0f + dot(u, u));
        const value_type invgam = 1.0f / gam;

        // unwind position one timestep backwards
        const Vec3 pos_0 = pos_1 - c * invgam * u;

        // --- masks (0.0f or 1.0f) ---
        // mask_skip: particle is in front of wall (no action needed)
        const value_type mask_skip = (pos_1(0) >= walloc) ? 1.0f : 0.0f;

        // close_enough: pre-push position was within c of the wall
        const value_type mask_close = (walloc0 - pos_0(0) <= c) ? 1.0f : 0.0f;

        // time fraction to wall crossing
        const value_type denom = betawall * c - c * u(0) * invgam;
        const value_type dt    = sstd::abs((pos_0(0) - walloc0) / (denom + EPS));

        // crossed: particle genuinely crossed the wall this timestep
        const value_type mask_crossed = (dt <= 1.0f) ? 1.0f : 0.0f;

        // composite masks (mutually exclusive, sum to 1)
        const value_type mask_refl = (1.0f - mask_skip) * mask_close * mask_crossed;
        const value_type mask_park = (1.0f - mask_skip) - mask_refl;

        // collision point (3D)
        const Vec3 pos_col = pos_0 + c * dt * invgam * u;

        // reflect x-velocity (general Lorentz boost; reduces to -ux when
        // betawall=0, gammawall=1)
        const value_type ux_new =
          gammawall * gammawall * gam *
          (2.0f * betawall - u(0) * invgam * (1.0f + betawall * betawall));

        const Vec3 u_new(ux_new, u(1), u(2));
        const value_type gam_new    = sstd::sqrt(1.0f + dot(u_new, u_new));
        const value_type invgam_new = 1.0f / gam_new;

        // reflected time fraction: remaining portion after collision
        const value_type dt_refl = sstd::min(
          1.0f,
          sstd::abs((pos_1(0) - pos_col(0)) / (pos_1(0) - pos_0(0) + EPS)));

        // new position after reflection
        const Vec3 pos_refl = pos_col + c * dt_refl * invgam_new * u_new;

        // --- safe current deposits with coordinate blending ---
        // For non-reflected particles, pos_col may be out-of-bounds.
        // Blend coordinates so non-reflected particles use p1l (always
        // valid in-tile), and zero the charge to produce no current.
        const Vec3 p1l = pos_1 - origo;

        // +q deposit: real current from pos_0 to pos_col
        const Vec3 dep_fwd_from = p1l + mask_refl * (pos_0 - pos_1);
        const Vec3 dep_fwd_to   = p1l + mask_refl * (pos_col - pos_1);
        zigzag_deposit_single(Jmds, dep_fwd_from, dep_fwd_to, mask_refl * charge);

        // -q deposit: cancel behind-wall current
        const Vec3 x1_deposit   = pos_refl - c * invgam_new * u_new;
        const Vec3 dep_rev_from = p1l + mask_refl * (x1_deposit - pos_1);
        const Vec3 dep_rev_to   = p1l + mask_refl * (pos_col - pos_1);
        zigzag_deposit_single(Jmds, dep_rev_from, dep_rev_to, mask_refl * (-charge));

        // --- final blending ---
        // x-position: 3-way blend
        pos_mds[idx][0] =
          mask_skip * pos_1(0) + mask_park * 1.0f + mask_refl * pos_refl(0);
        // y,z-position: only reflected particles change
        pos_mds[idx][1] = (1.0f - mask_refl) * pos_1(1) + mask_refl * pos_refl(1);
        pos_mds[idx][2] = (1.0f - mask_refl) * pos_1(2) + mask_refl * pos_refl(2);
        // x-velocity: skip keeps original, park zeroes, reflect sets ux_new
        vel_mds[idx][0] = mask_skip * u(0) + mask_refl * ux_new;
      })
    .wait();
}


/// Register a reflector wall with this tile.
template<std::size_t D>
void
  Tile<D>::register_reflector_wall(pic::reflector_wall wall)
{
  reflector_walls_.push_back(wall);
}

/// Reflect particles at all registered walls that overlap this tile.
///
/// Allocates a temporary correction-J grid, calls reflect_at_wall for every
/// particle species, then the correction is added to the tile's Yee lattice
/// current during the deposit phase.
template<std::size_t D>
void
  Tile<D>::reflect_particles()
{
  if(reflector_walls_.empty()) return;

  const auto wall_is_in_tile = [&](const reflector_wall& w) {
    return w.walloc >= value_type(this->mins[0]) - value_type(this->cfl_) &&
           w.walloc <= value_type(this->maxs[0]);
  };
  if(std::ranges::none_of(reflector_walls_, wall_is_in_tile)) return;

  // allocate and zero the correction J grid
  reflector_correction_J_.emplace(this->yee_lattice_.extents_with_halo());

  const auto genJmds = reflector_correction_J_->mds();
  tyvi::mdgrid_work {}
    .for_each_index(
      genJmds,
      [=](const auto idx, const auto tidx) { genJmds[idx][tidx] = 0; })
    .wait();

  const auto origo_pos = std::array { value_type(this->mins[0]) - this->halo_size,
                                      value_type(this->mins[1]) - this->halo_size,
                                      value_type(this->mins[2]) - this->halo_size };

  for(const auto& wall: reflector_walls_) {
    if(!wall_is_in_tile(wall)) continue;

    for(auto& [_, pbuff]: particle_buffs_) {
      pbuff
        .reflect_at_wall(wall, reflector_correction_J_.value(), origo_pos, this->cfl_);
    }
  }
}

/// Apply conducting-wall field BC: zero E_y and E_z behind each wall to
/// enforce the ideal-conductor constraint.
template<std::size_t D>
void
  Tile<D>::reflector_wall_field_bc()
{

  if(reflector_walls_.empty()) return;

  const auto nx = this->yee_lattice_.extents_with_halo()[0];

  for(const auto& wall: reflector_walls_) {
    if(wall.walloc >= value_type(this->maxs[0])) continue;

    // convert global wall location to lattice-local index (with halo offset)
    const auto wall_idx =
      static_cast<int>(wall.walloc - value_type(this->mins[0])) + this->halo_size;
    const auto wall_idx_clamped = static_cast<std::size_t>(std::max(wall_idx, 0));

    if(wall_idx_clamped < nx) {
      this->yee_lattice_.zero_transverse_E_behind_x(wall_idx_clamped);
    }
  }
}

/// Advance all wall positions by betawall * cfl each timestep.
template<std::size_t D>
void
  Tile<D>::advance_reflector_walls()
{
  if(reflector_walls_.empty()) return;

  for(auto& wall: reflector_walls_) {
    wall.walloc += wall.betawall * value_type(this->cfl_);
  }
}

}  // namespace pic


template class pic::Tile<3>;
