#include "core/emf/yee_lattice.h"
#include "core/pic/particle.h"
#include "core/pic/reflector_wall.h"
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

constexpr value_type EPS = 1e-10f;

}  // namespace


namespace pic {

/// Atomic zigzag current deposit for a single-particle sub-trajectory x1 → x2.
///
/// x1 and x2 are in lattice-local coordinates (relative to lattice origin).
/// Mirrors the atomic variant in particle_current_zigzag_1st.c++:276-329.
template<typename Jmds_t>
inline void zigzag_deposit_single(
  const Jmds_t Jmds,
  const Vec3 x1,
  const Vec3 x2,
  const value_type charge)
{
  const auto fi1 = Vec3(
    sstd::floor(x1(0)), sstd::floor(x1(1)), sstd::floor(x1(2)));
  const auto fi2 = Vec3(
    sstd::floor(x2(0)), sstd::floor(x2(1)), sstd::floor(x2(2)));

  const auto relay = [&](const runko::size_t j) -> value_type {
    const auto a  = sstd::min(fi1(j), fi2(j)) + value_type { 1 };
    const auto b1 = sstd::max(fi1(j), fi2(j));
    const auto b2 = value_type { 0.5 } * (x1(j) + x2(j));
    const auto b  = sstd::max(b1, b2);
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
    const auto si  = index.template as<runko::size_t>();
    auto* const Jx = &thrust::raw_reference_cast(Jmds[si.data][0]);
    auto* const Jy = &thrust::raw_reference_cast(Jmds[si.data][1]);
    auto* const Jz = &thrust::raw_reference_cast(Jmds[si.data][2]);

    sstd::atomic_add(Jx, current[0]);
    sstd::atomic_add(Jy, current[1]);
    sstd::atomic_add(Jz, current[2]);
  };

  // clang-format off
  store_current(i1,                   Vec3(Fx1*(1 - Wy1)*(1 - Wz1),
                                           Fy1*(1 - Wx1)*(1 - Wz1),
                                           Fz1*(1 - Wx1)*(1 - Wy1)));
  store_current(i2,                   Vec3(Fx2*(1 - Wy2)*(1 - Wz2),
                                           Fy2*(1 - Wx2)*(1 - Wz2),
                                           Fz2*(1 - Wx2)*(1 - Wy2)));
  store_current(i1 + Vec3uz(1, 0, 0), Vec3(0,
                                           Fy1*Wx1*(1 - Wz1),
                                           Fz1*Wx1*(1 - Wy1)));
  store_current(i2 + Vec3uz(1, 0, 0), Vec3(0,
                                           Fy2*Wx2*(1 - Wz2),
                                           Fz2*Wx2*(1 - Wy2)));
  store_current(i1 + Vec3uz(0, 1, 0), Vec3(Fx1*Wy1*(1 - Wz1),
                                           0,
                                           Fz1*(1 - Wx1)*Wy1));
  store_current(i2 + Vec3uz(0, 1, 0), Vec3(Fx2*Wy2*(1 - Wz2),
                                           0,
                                           Fz2*(1 - Wx2)*Wy2));
  store_current(i1 + Vec3uz(0, 0, 1), Vec3(Fx1*(1 - Wy1)*Wz1,
                                           Fy1*(1 - Wx1)*Wz1,
                                           0));
  store_current(i2 + Vec3uz(0, 0, 1), Vec3(Fx2*(1 - Wy2)*Wz2,
                                           Fy2*(1 - Wx2)*Wz2,
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
        // Position convention:
        //   pos_1    = current post-push position (may be behind wall)
        //   pos_0    = unwound pre-push position (start of timestep)
        //   pos_col  = wall collision point
        //   pos_refl = final reflected position (in front of wall)
        const Vec3 pos_1 = Vec3(pos_mds[idx]);

        // skip particles in front of the wall
        if(pos_1(0) >= walloc) return;

        const Vec3 u   = Vec3(vel_mds[idx]);
        const value_type gam = sstd::sqrt(1.0f + toolbox::dot(u, u));

        // unwind position one timestep backwards
        const Vec3 pos_0 = pos_1 - c*u/gam;

        // park artifact particle near domain start with zero x-velocity
        const auto park = [&] {
          pos_mds[idx][0] = 1.0f;
          vel_mds[idx][0] = 0;
        };

        // reject artifacts: particle was already far behind the wall
        if(walloc0 - pos_0(0) > c) { park(); return; }

        // time fraction to wall crossing
        const value_type denom = betawall * c - c * u(0)/gam;
        const value_type dt    = sstd::abs((pos_0(0) - walloc0) / (denom + EPS));

        // particle did not genuinely cross the wall this timestep
        if(dt > 1.0f) { park(); return; }

        // collision point (3D)
        const Vec3 pos_col = pos_0 + c*dt*u/gam;

        // dual deposit: +q from pos_0 to pos_col records real current before
        // reflection; -q from unwound-reflected-pos to pos_col (below) cancels
        // the spurious behind-wall current that normal deposit will add.
        zigzag_deposit_single(Jmds, pos_0 - origo, pos_col - origo, charge);

        // reflect x-velocity
        // stationary: ux' = -ux
        // moving:     ux' = gammawall^2 * gam * (2*beta - ux/gam * (1 + beta^2))
        const value_type ux_new =
          (betawall == 0 && gammawall == 1.0f)
            ? -u(0)
            : gammawall * gammawall * gam * (2.0f * betawall -
                 u(0) / gam * (1.0f + betawall * betawall));

        const Vec3 u_new(ux_new, u(1), u(2));
        const value_type gam_new    = sstd::sqrt(1.0f + toolbox::dot(u_new, u_new));
        const value_type invgam_new = 1.0f / gam_new;

        // reflected time fraction: remaining portion after collision
        const value_type dt_refl =
          sstd::min(1.0f,
                    sstd::abs((pos_1(0) - pos_col(0)) / (pos_1(0) - pos_0(0) + EPS)));

        // new position after reflection
        const Vec3 pos_refl = pos_col + c * dt_refl * invgam_new * u_new;

        // deposit -q current to cancel behind-wall current that normal deposit will add
        // normal deposit unwinds from pos_refl: x1_deposit = pos_refl - c * u_new / gam_new
        const Vec3 x1_deposit = pos_refl - c * invgam_new * u_new;
        zigzag_deposit_single(Jmds, x1_deposit - origo, pos_col - origo, -charge);

        // update all 3 dimensions: y/z also change because reflected trajectory
        // uses the new velocity over the remaining time fraction
        for(auto i = 0u; i < 3u; ++i) { pos_mds[idx][i] = pos_refl(i); }
        vel_mds[idx][0] = ux_new;
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
    return w.walloc >= value_type(this->mins[0]) - value_type(this->cfl_)
        && w.walloc <= value_type(this->maxs[0]);
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

  const auto origo_pos =
    std::array { value_type(this->mins[0]) - this->halo_size,
                 value_type(this->mins[1]) - this->halo_size,
                 value_type(this->mins[2]) - this->halo_size };

  for(const auto& wall: reflector_walls_) {
    if(!wall_is_in_tile(wall)) continue;

    for(auto& [_, pbuff]: particle_buffs_) {
      pbuff.reflect_at_wall(wall, reflector_correction_J_.value(), origo_pos, this->cfl_);
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
      static_cast<int>(wall.walloc - value_type(this->mins[0])) +
      this->halo_size;
    const auto wall_idx_clamped =
      static_cast<runko::size_t>(std::max(wall_idx, 0));

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
