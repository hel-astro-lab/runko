#pragma once

#include "tyvi/mdgrid.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <concepts>
#include <cstddef>
#include <experimental/mdspan>
#include <tuple>
#include <type_traits>
#include <vector>

namespace emf2 {

struct YeeLatticeCtorArgs {
  std::size_t halo_size {};
  std::size_t Nx {}, Ny {}, Nz {};
};

struct YeeLatticeFieldsAtPoint {
  double Ex, Ey, Ez;
  double Bx, By, Bz;
  double Jx, Jy, Jz;
};

/// Concept for functions that return Yee Lattice fields at some grid point.
///
/// F has sematic requirement that fields have to be calculated at correct positions
/// in the Yee Lattice.
template<typename F>
concept yee_lattice_fields_function =
  std::invocable<F, std::size_t, std::size_t, std::size_t> and
  std::same_as<
    std::invoke_result_t<F, std::size_t, std::size_t, std::size_t>,
    YeeLatticeFieldsAtPoint>;

/// Yee lattice of plasma quantities in tyvi::mdgrid continers.
class [[nodiscard]] YeeLattice {
public:
  using grid_extents_type = std::dextents<std::size_t, 3>;
  static constexpr auto vec_element =
    tyvi::mdgrid_element_descriptor<float> { .rank = 1, .dim = 3 };

  using vec_grid = tyvi::mdgrid<vec_element, grid_extents_type>;

  using YeeLatticeHostCopy = tyvi::mdgrid_buffer<
    std::vector<YeeLatticeFieldsAtPoint>,
    std::extents<std::size_t>,
    std::layout_right,
    vec_grid::grid_extents_type,
    vec_grid::grid_layout_type>;

private:
  std::size_t halo_size_;
  std::array<std::size_t, 3> extents_wout_halo_;


  /// Electric field
  vec_grid E_;

  /// Magnetic field
  vec_grid B_;

  /// Current
  vec_grid J_;

public:
  explicit YeeLattice(YeeLatticeCtorArgs);

  /* mds getters are implemented in header, due to not wanting to write return type. */

  /// Returns tuple of mdspans to yee_lattice_ staging buffer {E, B, J} with halo.
  [[nodiscard]] auto staging_mds_with_halo()
  {
    return std::tuple { E_.staging_mds(), B_.staging_mds(), J_.staging_mds() };
  }

  /// Returns tuple of mdspans to yee_lattice_ staging buffer {E, B, J} without halo.
  [[nodiscard]] auto staging_mds_wout_halo()
  {
    auto [Emds, Bmds, Jmds] = staging_mds_with_halo();

    const auto x = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[0] };
    const auto y = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[1] };
    const auto z = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[2] };

    return std::tuple { std::submdspan(std::move(Emds), x, y, z),
                        std::submdspan(std::move(Bmds), x, y, z),
                        std::submdspan(std::move(Jmds), x, y, z) };
  }

  /// Returns tuple of mdspans to yee_lattice_ device buffer {E, B, J} with halo.
  [[nodiscard]] auto mds_with_halo()
  {
    return std::tuple { E_.mds(), B_.mds(), J_.mds() };
  }

  /// Returns tuple of mdspans to yee_lattice_ device buffer {E, B, J} without halo.
  [[nodiscard]] auto mds_wout_halo()
  {
    auto [Emds, Bmds, Jmds] = mds_with_halo();

    const auto x = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[0] };
    const auto y = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[1] };
    const auto z = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[2] };

    return std::tuple { std::submdspan(std::move(Emds), x, y, z),
                        std::submdspan(std::move(Bmds), x, y, z),
                        std::submdspan(std::move(Jmds), x, y, z) };
  }


  [[nodiscard]] std::array<std::size_t, 3> extents_wout_halo() const;
  [[nodiscard]] std::array<std::size_t, 3> extents_with_halo() const;
  [[nodiscard]] std::size_t halo_size() const;


  /// Initializes E, B and J in non-halo region.
  void set_EBJ(yee_lattice_fields_function auto&& f)
  {
    const auto [Emds, Bmds, Jmds] = staging_mds_wout_halo();

    for(const auto idx: tyvi::sstd::index_space(Emds)) {
      const auto F = f(idx[0], idx[1], idx[2]);
      Emds[idx][0] = F.Ex;
      Emds[idx][1] = F.Ey;
      Emds[idx][2] = F.Ez;
      Bmds[idx][0] = F.Bx;
      Bmds[idx][1] = F.By;
      Bmds[idx][2] = F.Bz;
      Jmds[idx][0] = F.Jx;
      Jmds[idx][1] = F.Jy;
      Jmds[idx][2] = F.Jz;
    }

    auto wE = tyvi::mdgrid_work {};
    auto wB = tyvi::mdgrid_work {};
    auto wJ = tyvi::mdgrid_work {};

    auto wE1 = wE.sync_from_staging(E_);
    auto wB1 = wB.sync_from_staging(B_);
    auto wJ1 = wJ.sync_from_staging(J_);

    tyvi::when_all(wE1, wB1, wJ1).wait();
  }

  /// Get copy of E, B and J in non-halo region.
  YeeLatticeHostCopy get_EBJ();
};

}  // namespace emf2
