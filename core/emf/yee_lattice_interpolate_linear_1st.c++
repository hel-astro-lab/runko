#include "core/emf/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include <array>
#include <cstddef>

namespace {

/// Trilinear interpolation of 8 corner values at fractional position (dx, dy, dz).
/// Arguments v_{ic,jc,kc} for ic,jc,kc in {0,1}, ordered as
///   v000, v100, v010, v110, v001, v101, v011, v111.
constexpr auto
  lerp3D(const auto dx, const auto dy, const auto dz,
         const auto v000, const auto v100, const auto v010, const auto v110,
         const auto v001, const auto v101, const auto v011, const auto v111)
{
  const auto c00 = v000 + dx * (v100 - v000);
  const auto c10 = v010 + dx * (v110 - v010);
  const auto c01 = v001 + dx * (v101 - v001);
  const auto c11 = v011 + dx * (v111 - v011);

  const auto c0 = c00 + dy * (c10 - c00);
  const auto c1 = c01 + dy * (c11 - c01);

  return c0 + dz * (c1 - c0);
}

}  // namespace

namespace emf {

YeeLattice::InterpolatedEB
  YeeLattice::interpolate_EB_linear_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const runko::VecList<value_type>& pos) const
{
  tyvi::mdgrid_work w {};
  return interpolate_EB_linear_1st(w, lattice_origo_coordinates, pos);
}

YeeLattice::InterpolatedEB
  YeeLattice::interpolate_EB_linear_1st(
    const tyvi::mdgrid_work& w,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const runko::VecList<value_type>& pos) const
{
  auto Eipo = runko::VecList<value_type>(pos.extents());
  auto Bipo = runko::VecList<value_type>(pos.extents());

  const auto posmds   = pos.mds();
  const auto Eipo_mds = Eipo.mds();
  const auto Bipo_mds = Bipo.mds();

  // Note that these include the halo.
  const auto Emds = E_.mds();
  const auto Bmds = B_.mds();

  // Fill in Eipo_mds and Bipo_mds via unrolled 8-point gather + trilinear interp.
  w.for_each_index(
     posmds,
     [=](const auto idx) {
       const auto pos_vec        = toolbox::Vec3<value_type>(posmds[idx]);
       const auto origo          = toolbox::Vec3<value_type>(lattice_origo_coordinates);
       const auto pos_in_lattice = pos_vec - origo;

       // This should be equivelant to floor, as we assume that the particles
       // are inside the lattice, such that the subtractions are always positive.
       const auto index_in_lattice = pos_in_lattice.template as<std::size_t>();
       const auto [i, j, k]       = index_in_lattice.data;

       const auto dp = pos_in_lattice - index_in_lattice.template as<value_type>();
       const auto dx = dp[0], dy = dp[1], dz = dp[2];

       static constexpr auto mean2 = [](const auto a, const auto b) {
         return (a + b) * value_type { 0.5 };
       };
       static constexpr auto mean4 =
         [](const auto a, const auto b, const auto c, const auto d) {
           return (a + b + c + d) * value_type { 0.25 };
         };

       // Ex: face-centered along x, average adjacent x-cells
       Eipo_mds[idx][0] = lerp3D(dx, dy, dz,
         mean2(Emds[i-1, j  , k  ][0], Emds[i  , j  , k  ][0]),
         mean2(Emds[i  , j  , k  ][0], Emds[i+1, j  , k  ][0]),
         mean2(Emds[i-1, j+1, k  ][0], Emds[i  , j+1, k  ][0]),
         mean2(Emds[i  , j+1, k  ][0], Emds[i+1, j+1, k  ][0]),
         mean2(Emds[i-1, j  , k+1][0], Emds[i  , j  , k+1][0]),
         mean2(Emds[i  , j  , k+1][0], Emds[i+1, j  , k+1][0]),
         mean2(Emds[i-1, j+1, k+1][0], Emds[i  , j+1, k+1][0]),
         mean2(Emds[i  , j+1, k+1][0], Emds[i+1, j+1, k+1][0]));

       // Ey: face-centered along y, average adjacent y-cells
       Eipo_mds[idx][1] = lerp3D(dx, dy, dz,
         mean2(Emds[i  , j-1, k  ][1], Emds[i  , j  , k  ][1]),
         mean2(Emds[i+1, j-1, k  ][1], Emds[i+1, j  , k  ][1]),
         mean2(Emds[i  , j  , k  ][1], Emds[i  , j+1, k  ][1]),
         mean2(Emds[i+1, j  , k  ][1], Emds[i+1, j+1, k  ][1]),
         mean2(Emds[i  , j-1, k+1][1], Emds[i  , j  , k+1][1]),
         mean2(Emds[i+1, j-1, k+1][1], Emds[i+1, j  , k+1][1]),
         mean2(Emds[i  , j  , k+1][1], Emds[i  , j+1, k+1][1]),
         mean2(Emds[i+1, j  , k+1][1], Emds[i+1, j+1, k+1][1]));

       // Ez: face-centered along z, average adjacent z-cells
       Eipo_mds[idx][2] = lerp3D(dx, dy, dz,
         mean2(Emds[i  , j  , k-1][2], Emds[i  , j  , k  ][2]),
         mean2(Emds[i+1, j  , k-1][2], Emds[i+1, j  , k  ][2]),
         mean2(Emds[i  , j+1, k-1][2], Emds[i  , j+1, k  ][2]),
         mean2(Emds[i+1, j+1, k-1][2], Emds[i+1, j+1, k  ][2]),
         mean2(Emds[i  , j  , k  ][2], Emds[i  , j  , k+1][2]),
         mean2(Emds[i+1, j  , k  ][2], Emds[i+1, j  , k+1][2]),
         mean2(Emds[i  , j+1, k  ][2], Emds[i  , j+1, k+1][2]),
         mean2(Emds[i+1, j+1, k  ][2], Emds[i+1, j+1, k+1][2]));

       // Bx: edge-centered along x, average 4 adjacent cells in y-z plane
       Bipo_mds[idx][0] = lerp3D(dx, dy, dz,
         mean4(Bmds[i  , j  , k  ][0], Bmds[i  , j-1, k  ][0],
               Bmds[i  , j  , k-1][0], Bmds[i  , j-1, k-1][0]),
         mean4(Bmds[i+1, j  , k  ][0], Bmds[i+1, j-1, k  ][0],
               Bmds[i+1, j  , k-1][0], Bmds[i+1, j-1, k-1][0]),
         mean4(Bmds[i  , j+1, k  ][0], Bmds[i  , j  , k  ][0],
               Bmds[i  , j+1, k-1][0], Bmds[i  , j  , k-1][0]),
         mean4(Bmds[i+1, j+1, k  ][0], Bmds[i+1, j  , k  ][0],
               Bmds[i+1, j+1, k-1][0], Bmds[i+1, j  , k-1][0]),
         mean4(Bmds[i  , j  , k+1][0], Bmds[i  , j-1, k+1][0],
               Bmds[i  , j  , k  ][0], Bmds[i  , j-1, k  ][0]),
         mean4(Bmds[i+1, j  , k+1][0], Bmds[i+1, j-1, k+1][0],
               Bmds[i+1, j  , k  ][0], Bmds[i+1, j-1, k  ][0]),
         mean4(Bmds[i  , j+1, k+1][0], Bmds[i  , j  , k+1][0],
               Bmds[i  , j+1, k  ][0], Bmds[i  , j  , k  ][0]),
         mean4(Bmds[i+1, j+1, k+1][0], Bmds[i+1, j  , k+1][0],
               Bmds[i+1, j+1, k  ][0], Bmds[i+1, j  , k  ][0]));

       // By: edge-centered along y, average 4 adjacent cells in x-z plane
       Bipo_mds[idx][1] = lerp3D(dx, dy, dz,
         mean4(Bmds[i  , j  , k  ][1], Bmds[i-1, j  , k  ][1],
               Bmds[i  , j  , k-1][1], Bmds[i-1, j  , k-1][1]),
         mean4(Bmds[i+1, j  , k  ][1], Bmds[i  , j  , k  ][1],
               Bmds[i+1, j  , k-1][1], Bmds[i  , j  , k-1][1]),
         mean4(Bmds[i  , j+1, k  ][1], Bmds[i-1, j+1, k  ][1],
               Bmds[i  , j+1, k-1][1], Bmds[i-1, j+1, k-1][1]),
         mean4(Bmds[i+1, j+1, k  ][1], Bmds[i  , j+1, k  ][1],
               Bmds[i+1, j+1, k-1][1], Bmds[i  , j+1, k-1][1]),
         mean4(Bmds[i  , j  , k+1][1], Bmds[i-1, j  , k+1][1],
               Bmds[i  , j  , k  ][1], Bmds[i-1, j  , k  ][1]),
         mean4(Bmds[i+1, j  , k+1][1], Bmds[i  , j  , k+1][1],
               Bmds[i+1, j  , k  ][1], Bmds[i  , j  , k  ][1]),
         mean4(Bmds[i  , j+1, k+1][1], Bmds[i-1, j+1, k+1][1],
               Bmds[i  , j+1, k  ][1], Bmds[i-1, j+1, k  ][1]),
         mean4(Bmds[i+1, j+1, k+1][1], Bmds[i  , j+1, k+1][1],
               Bmds[i+1, j+1, k  ][1], Bmds[i  , j+1, k  ][1]));

       // Bz: edge-centered along z, average 4 adjacent cells in x-y plane
       Bipo_mds[idx][2] = lerp3D(dx, dy, dz,
         mean4(Bmds[i  , j  , k  ][2], Bmds[i-1, j  , k  ][2],
               Bmds[i  , j-1, k  ][2], Bmds[i-1, j-1, k  ][2]),
         mean4(Bmds[i+1, j  , k  ][2], Bmds[i  , j  , k  ][2],
               Bmds[i+1, j-1, k  ][2], Bmds[i  , j-1, k  ][2]),
         mean4(Bmds[i  , j+1, k  ][2], Bmds[i-1, j+1, k  ][2],
               Bmds[i  , j  , k  ][2], Bmds[i-1, j  , k  ][2]),
         mean4(Bmds[i+1, j+1, k  ][2], Bmds[i  , j+1, k  ][2],
               Bmds[i+1, j  , k  ][2], Bmds[i  , j  , k  ][2]),
         mean4(Bmds[i  , j  , k+1][2], Bmds[i-1, j  , k+1][2],
               Bmds[i  , j-1, k+1][2], Bmds[i-1, j-1, k+1][2]),
         mean4(Bmds[i+1, j  , k+1][2], Bmds[i  , j  , k+1][2],
               Bmds[i+1, j-1, k+1][2], Bmds[i  , j-1, k+1][2]),
         mean4(Bmds[i  , j+1, k+1][2], Bmds[i-1, j+1, k+1][2],
               Bmds[i  , j  , k+1][2], Bmds[i-1, j  , k+1][2]),
         mean4(Bmds[i+1, j+1, k+1][2], Bmds[i  , j+1, k+1][2],
               Bmds[i+1, j  , k+1][2], Bmds[i  , j  , k+1][2]));
     })
    .wait();
  return { std::move(Eipo), std::move(Bipo) };
}

}  // namespace emf
