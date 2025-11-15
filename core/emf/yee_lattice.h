#pragma once

#include "core/mdgrid_common.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <experimental/mdspan>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace emf {

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
  using value_type = float;
  using VecGrid    = runko::VecGrid<value_type>;
  using VecGridMDS = std::remove_cvref_t<decltype(std::declval<VecGrid&>().mds())>;

  using YeeLatticeHostCopy = tyvi::mdgrid_buffer<
    std::vector<YeeLatticeFieldsAtPoint>,
    std::extents<std::size_t>,
    std::layout_right,
    VecGrid::grid_extents_type,
    VecGrid::grid_layout_type>;


  /// Lattice consists of  27 different regions which are labeled with {i, j, k}.
  ///
  /// Valid dir labels are -1, 0 and 1 and specify coordinate extents for one dimension.
  /// These correspond to regions [0, halo_size_), [halo_size_, halo_size + extent)
  /// and [halo_size + extent, 2 * halo_size + extent) respectively.
  using dir_type = std::array<int, 3>;

  [[nodiscard]] static constexpr dir_type invert_dir(const dir_type& dir)
  {
    return dir_type { -dir[0], -dir[1], -dir[2] };
  }

private:
  std::size_t halo_size_;
  std::array<std::size_t, 3> extents_wout_halo_;


  /// Electric field
  VecGrid E_;

  /// Magnetic field
  VecGrid B_;

  /// Current
  VecGrid J_;

  /* FIXME: Use std::integral_constant when possible in mds helpers below. */

  /// Convinience function to get non-halo region of some field.
  ///
  /// This could have variadic number of arguments and be used as:
  ///
  /// const auto [Emds, Bmds, Jmds] = nonhalo_submds(E_.mds(), B_.mds(), J_.mds());
  ///
  /// But Cray Clang on Hile (CC) with OpenMP does not support capturing
  /// structured bindings to lambdas, so above has to be written as:
  ///
  /// const auto Emds = nonhalo_submds(E_.mds());
  /// const auto Bmds = nonhalo_submds(B_.mds());
  /// const auto Jmds = nonhalo_submds(J_.mds());
  template<typename MDS>
  [[nodiscard]]
  auto nonhalo_submds(MDS&& mds);

  /// Throws if given mdspan's extents differ from extents_with_halo.
  void assert_mds_spans_whole_lattice(const auto& mds) const;

  /// Get subregion from whole field mds specified by dir.
  ///
  /// mds is assumed to span the whole lattice including non-halo regions.
  template<typename MDS>
  [[nodiscard]]
  auto subregion(dir_type dir, MDS&& mds) const;

  /// Get region in non-halo region corresponding to halo region from `subregion`.
  ///
  /// mds is assumed to span the whole lattice including non-halo regions.
  template<typename MDS>
  [[nodiscard]]
  auto corresponding_subregion(dir_type dir, MDS&& mds) const;

public:
  explicit YeeLattice(YeeLatticeCtorArgs);

  [[nodiscard]] std::array<std::size_t, 3> extents_wout_halo() const;
  [[nodiscard]] std::array<std::size_t, 3> extents_with_halo() const;
  [[nodiscard]] std::size_t halo_size() const;

  /// Initializes E, B and J in non-halo region.
  void set_EBJ(yee_lattice_fields_function auto&& f);

  /// Get copy of E, B and J in non-halo region.
  YeeLatticeHostCopy get_EBJ();

  /// Get copy of E, B and J including halo regions.
  YeeLatticeHostCopy get_EBJ_with_halo();

  /// Advance B by half time step using FDTD2 scheme in non-halo region.
  void push_b_FDTD2(value_type dt);

  /// Advance E by full time step using FDTD2 scheme in non-halo region.
  void push_e_FDTD2(value_type dt);

  /// E -= J in non-halo region.
  void subtract_J_from_E();

  /// Advance B by half time step using FDTD2 scheme in non-halo region asynchronously.
  void push_b_FDTD2(const tyvi::mdgrid_work&, value_type dt);

  /// Advance E by full time step using FDTD2 scheme in non-halo region asynchronously.
  void push_e_FDTD2(const tyvi::mdgrid_work&, value_type dt);

  /// E -= J in non-halo region asynchronously.
  void subtract_J_from_E(const tyvi::mdgrid_work&);

  [[nodiscard]]
  auto span_E() &;

  [[nodiscard]]
  auto span_E() const&;

  [[nodiscard]]
  auto span_B() &;

  [[nodiscard]]
  auto span_B() const&;

  [[nodiscard]]
  auto span_J() &;

  [[nodiscard]]
  auto span_J() const&;


  /// Set E subregion specified by dir from other tile.
  void set_E_in_subregion(dir_type dir, const YeeLattice& other);

  /// Set B subregion specified by dir from other tile.
  void set_B_in_subregion(dir_type dir, const YeeLattice& other);

  /// Set J subregion specified by dir from other tile.
  void set_J_in_subregion(dir_type dir, const YeeLattice& other);

  /// Set E subregion specified by dir from other tile asynchronously.
  void
    set_E_in_subregion(const tyvi::mdgrid_work&, dir_type dir, const YeeLattice& other);

  /// Set B subregion specified by dir from other tile asynchronously.
  void
    set_B_in_subregion(const tyvi::mdgrid_work&, dir_type dir, const YeeLattice& other);

  /// Set J subregion specified by dir from other tile asynchronously.
  void
    set_J_in_subregion(const tyvi::mdgrid_work&, dir_type dir, const YeeLattice& other);

  /// Add J from other's subregion to corresponding region in this.
  ///
  /// Note that compared to set_{E,B,J}_in_subregion,
  /// this reads other's halo and adds it to non-halo of this,
  /// while the other's read non-halo and add it to halo of this.
  void add_to_J_from_subregion(dir_type dir, const YeeLattice& other);

  /// Add J from other's subregion to corresponding region in this asynchronously.
  ///
  /// Note that compared to set_{E,B,J}_in_subregion,
  /// this reads other's halo and adds it to non-halo of this,
  /// while the other's read non-halo and add it to halo of this.
  void add_to_J_from_subregion(
    const tyvi::mdgrid_work&,
    dir_type dir,
    const YeeLattice& other);

  struct [[nodiscard]] InterpolatedEB {
    runko::VecList<value_type> E, B;
  };

  /// Interpolate E and B to given coordinates using linear_1st interpolation.
  ///
  /// Lattice index (0, 0, 0) is interperted to be at lattice_origo_coordinates.
  /// Note that (0, 0, 0) is in the halo region.
  InterpolatedEB interpolate_EB_linear_1st(
    std::array<value_type, 3> lattice_origo_coordinates,
    const runko::VecList<value_type>& coordinates) const;

  /// Interpolate E and B to given coordinates using linear_1st interpolation (async).
  ///
  /// Lattice index (0, 0, 0) is interperted to be at lattice_origo_coordinates.
  /// Note that (0, 0, 0) is in the halo region.
  ///
  /// By async means that work is executed on the given work which is synchronized
  /// before returning.
  InterpolatedEB interpolate_EB_linear_1st(
    const tyvi::mdgrid_work&,
    std::array<value_type, 3> lattice_origo_coordinates,
    const runko::VecList<value_type>& coordinates) const;

  /// Set J = 0 everywhere including in halo.
  void clear_current();

  /// Set J = 0 everywhere including in halo asynchronously.
  void clear_current(const tyvi::mdgrid_work&);

  /// Represents a set of locations and corresponding currents.
  struct [[nodiscard]] CurrentContributions {
    thrust::device_vector<std::array<std::size_t, 3>> locations;
    thrust::device_vector<std::array<value_type, 3>> currents;
  };

  /// Add given currents to J.
  ///
  /// It is assumed that every location appears at most once and
  /// that the location indices include the halo regions.
  void deposit_current(const CurrentContributions&);

  /// Add given currents to J.
  ///
  /// Throws if given grid is not same size as the lattice with halo.
  void deposit_current(const runko::VecGrid<value_type>&);

  /// Add given currents to J asynchronously.
  ///
  /// It is assumed that every location appears at most once and
  /// that the location indices include the halo regions.
  void deposit_current(const tyvi::mdgrid_work&, const CurrentContributions&);

  /// Add given currents to J asynchronously.
  ///
  /// Throws if given grid is not same size as the lattice with halo.
  void deposit_current(const tyvi::mdgrid_work&, const runko::VecGrid<value_type>&);

  /// Apply digital 2nd order one-pass binomial filter for J.
  void filter_current_binomial2();

  /// Apply digital 2nd order one-pass binomial filter for J asynchronouosly.
  ///
  /// Asynchronously  means that work is executed on the given work which is
  /// synchronized before returning.
  void filter_current_binomial2(const tyvi::mdgrid_work&);

  /// Returns mdspans to host accessible E, B and J in non-halo region.
  auto view_EBJ_on_host();

  [[nodiscard]] YeeLattice::VecGridMDS::mapping_type
    grid_mapping_with_halo() const noexcept;
};

inline auto
  YeeLattice::span_E() &
{
  return this->E_.span();
}

inline auto
  YeeLattice::span_E() const&
{
  return this->E_.span();
}

inline auto
  YeeLattice::span_B() &
{
  return this->B_.span();
}

inline auto
  YeeLattice::span_B() const&
{
  return this->B_.span();
}

inline auto
  YeeLattice::span_J() &
{
  return this->J_.span();
}

inline auto
  YeeLattice::span_J() const&
{
  return this->J_.span();
}

template<typename MDS>
auto
  YeeLattice::nonhalo_submds(MDS&& mds)
{
  const auto x = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[0] };
  const auto y = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[1] };
  const auto z = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[2] };

  return std::submdspan(std::forward<MDS>(mds), x, y, z);
}

void
  YeeLattice::assert_mds_spans_whole_lattice(const auto& mds) const
{
  const auto mds_extents =
    std::views::iota(0uz, mds.rank()) |
    std::views::transform([&](const std::size_t i) { return mds.extent(i); });


  if(not std::ranges::equal(mds_extents, extents_with_halo())) {
    throw std::runtime_error { "Given mdspan's extents != extents_with_halo." };
  }
}

template<typename MDS>
auto
  YeeLattice::subregion(dir_type dir, MDS&& mds) const
{
  assert_mds_spans_whole_lattice(mds);

  auto oneD_dir_to_index_extent =
    [&, this](const std::size_t i) -> std::tuple<std::size_t, std::size_t> {
    switch(dir[i]) {
      case -1: return { 0uz, halo_size_ };
      case 0: return { halo_size_, halo_size_ + extents_wout_halo_[i] };
      case 1:
        return { halo_size_ + extents_wout_halo_[i],
                 2uz * halo_size_ + extents_wout_halo_[i] };
      default:
        throw std::logic_error { std::format("dir[{}] = {} != -1, 0 or 1", i, dir[i]) };
    }
  };

  const auto x = oneD_dir_to_index_extent(0);
  const auto y = oneD_dir_to_index_extent(1);
  const auto z = oneD_dir_to_index_extent(2);

  return std::submdspan(std::forward<MDS>(mds), x, y, z);
}

inline void
  YeeLattice::set_EBJ(yee_lattice_fields_function auto&& f)
{
  const auto [Emds, Bmds, Jmds] = std::tuple { nonhalo_submds(E_.staging_mds()),
                                               nonhalo_submds(B_.staging_mds()),
                                               nonhalo_submds(J_.staging_mds()) };

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

  wE.sync_from_staging(E_);
  wB.sync_from_staging(B_);
  wJ.sync_from_staging(J_);

  tyvi::when_all(wE, wB, wJ);
  wE.wait();
}

inline auto
  YeeLattice::view_EBJ_on_host()
{
  auto wE = tyvi::mdgrid_work {};
  auto wB = tyvi::mdgrid_work {};
  auto wJ = tyvi::mdgrid_work {};

  wE.sync_to_staging(E_);
  wB.sync_to_staging(B_);
  wJ.sync_to_staging(J_);

  tyvi::when_all(wE, wB, wJ);
  wE.wait();

  const auto Emds = nonhalo_submds(E_.staging_mds());
  const auto Bmds = nonhalo_submds(B_.staging_mds());
  const auto Jmds = nonhalo_submds(J_.staging_mds());

  return std::tuple { Emds, Bmds, Jmds };
}

template<typename MDS>
inline auto
  YeeLattice::corresponding_subregion(dir_type dir, MDS&& mds) const
{
  assert_mds_spans_whole_lattice(mds);

  auto oneD_dir_to_index_extent =
    [&, this](const std::size_t i) -> std::tuple<std::size_t, std::size_t> {
    switch(invert_dir(dir)[i]) {
      case -1: return { halo_size_, 2u * halo_size_ };
      case 0: return { halo_size_, halo_size_ + extents_wout_halo_[i] };
      case 1: return { extents_wout_halo_[i], halo_size_ + extents_wout_halo_[i] };
      default:
        throw std::logic_error { std::format("dir[{}] = {} != -1, 0 or 1", i, dir[i]) };
    }
  };

  const auto x = oneD_dir_to_index_extent(0);
  const auto y = oneD_dir_to_index_extent(1);
  const auto z = oneD_dir_to_index_extent(2);

  return std::submdspan(std::forward<MDS>(mds), x, y, z);
}

}  // namespace emf
