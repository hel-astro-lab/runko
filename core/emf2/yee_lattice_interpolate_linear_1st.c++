#include "core/emf2/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace {

using lerp3D_extents = std::extents<std::size_t, 2, 2, 2>;

/// Linear interpolation of mds values to pos.
///
/// Values in mds are assumed to be located at (i, j, k),
/// where i, j and k are either 0 or 1. Pos is assumed to be
/// inside the unit cube and measured from (0, 0, 0).
constexpr auto
  lerp3D(const auto& pos, const auto& mds)
{
  using MDS = std::remove_cvref_t<decltype(mds)>;
  static_assert(MDS::rank() == 3);
  static_assert(MDS::static_extent(0) == 2);
  static_assert(MDS::static_extent(1) == 2);
  static_assert(MDS::static_extent(2) == 2);

  constexpr auto lerp = [](const auto x, const auto A, const auto B) {
    return (1 - x) * A + x * B;
  };

  const auto c00 = lerp(pos[0], mds[0, 0, 0], mds[1, 0, 0]);
  const auto c01 = lerp(pos[0], mds[0, 0, 1], mds[1, 0, 1]);
  const auto c10 = lerp(pos[0], mds[0, 1, 0], mds[1, 1, 0]);
  const auto c11 = lerp(pos[0], mds[0, 1, 1], mds[1, 1, 1]);

  const auto c0 = lerp(pos[1], c00, c10);
  const auto c1 = lerp(pos[1], c01, c11);

  return lerp(pos[2], c0, c1);
}


template<std::signed_integral To>
constexpr To
  floor(const std::floating_point auto x)
{
  const auto y      = static_cast<To>(x);
  const auto is_neg = static_cast<To>(x < 0);
  return y - is_neg;
};
}  // namespace

namespace emf2 {

YeeLattice::InterpolatedEB
  YeeLattice::interpolate_EB_linear_1st(
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

  // Fill in Eipo_mds.
  tyvi::mdgrid_work {}
    .for_each_index(
      posmds,
      [=](const auto idx) {
        const auto i =
          floor<std::size_t>(posmds[idx][0] - lattice_origo_coordinates[0]);
        const auto j =
          floor<std::size_t>(posmds[idx][1] - lattice_origo_coordinates[1]);
        const auto k =
          floor<std::size_t>(posmds[idx][2] - lattice_origo_coordinates[2]);

        auto Ex_means = std::array<value_type, 8> {};
        auto Ey_means = std::array<value_type, 8> {};
        auto Ez_means = std::array<value_type, 8> {};

        auto Ex_means_mds = std::mdspan<value_type, lerp3D_extents>(Ex_means.data());
        auto Ey_means_mds = std::mdspan<value_type, lerp3D_extents>(Ey_means.data());
        auto Ez_means_mds = std::mdspan<value_type, lerp3D_extents>(Ez_means.data());

        auto Bx_means = std::array<value_type, 8> {};
        auto By_means = std::array<value_type, 8> {};
        auto Bz_means = std::array<value_type, 8> {};

        auto Bx_means_mds = std::mdspan<value_type, lerp3D_extents>(Bx_means.data());
        auto By_means_mds = std::mdspan<value_type, lerp3D_extents>(By_means.data());
        auto Bz_means_mds = std::mdspan<value_type, lerp3D_extents>(Bz_means.data());

        static constexpr auto mean = []<typename... T>(const T... vals) {
          const auto sum = (vals + ...);
          return sum / static_cast<std::remove_cvref_t<decltype(sum)>>(sizeof...(vals));
        };

        for(const auto [ic, jc, kc]: tyvi::sstd::index_space(Ex_means_mds)) {
          const auto [ii, jj, kk] =
            std::array<std::size_t, 3> { static_cast<std::size_t>(i + ic),
                                         static_cast<std::size_t>(j + jc),
                                         static_cast<std::size_t>(k + kc) };

          Ex_means_mds[ic, jc, kc] = mean(Emds[ii - 1, jj, kk][0], Emds[ii, jj, kk][0]);
          Ey_means_mds[ic, jc, kc] = mean(Emds[ii, jj - 1, kk][1], Emds[ii, jj, kk][1]);
          Ez_means_mds[ic, jc, kc] = mean(Emds[ii, jj, kk - 1][2], Emds[ii, jj, kk][2]);

          Bx_means_mds[ic, jc, kc] = mean(
            Bmds[ii, jj, kk][0],
            Bmds[ii, jj - 1, kk][0],
            Bmds[ii, jj, kk - 1][0],
            Bmds[ii, jj - 1, kk - 1][0]);

          By_means_mds[ic, jc, kc] = mean(
            Bmds[ii, jj, kk][1],
            Bmds[ii - 1, jj, kk][1],
            Bmds[ii, jj, kk - 1][1],
            Bmds[ii - 1, jj, kk - 1][1]);

          Bz_means_mds[ic, jc, kc] = mean(
            Bmds[ii, jj, kk][2],
            Bmds[ii - 1, jj, kk][2],
            Bmds[ii, jj - 1, kk][2],
            Bmds[ii - 1, jj - 1, kk][2]);
        }

        const auto delta_pos =
          std::array<value_type, 3> { static_cast<value_type>(posmds[idx][0] - i),
                                      static_cast<value_type>(posmds[idx][1] - j),
                                      static_cast<value_type>(posmds[idx][2] - k) };

        Eipo_mds[idx][0] = lerp3D(delta_pos, Ex_means_mds);
        Eipo_mds[idx][1] = lerp3D(delta_pos, Ey_means_mds);
        Eipo_mds[idx][2] = lerp3D(delta_pos, Ez_means_mds);

        Bipo_mds[idx][0] = lerp3D(delta_pos, Bx_means_mds);
        Bipo_mds[idx][1] = lerp3D(delta_pos, By_means_mds);
        Bipo_mds[idx][2] = lerp3D(delta_pos, Bz_means_mds);
      })
    .wait();
  return { std::move(Eipo), std::move(Bipo) };
}

}  // namespace emf2
