#include "core/emf/yee_lattice.h"

#include "thrust/execution_policy.h"
#include "thrust/for_each.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/reduce.h"
#include "thrust/tuple.h"
#include "tyvi/mdspan.h"

#include <algorithm>
#include <numbers>
#include <ranges>
#include <sstream>
#include <string>

namespace {

template<typename MDSfrom, typename MDSto>
void
  mds_copy(const tyvi::mdgrid_work& w, MDSfrom&& from, MDSto&& to)
{

  if(from.extents() != to.extents()) {
    throw std::runtime_error { "Can not copy from/to different sized mdspan." };
  }

  w.for_each_index(from, [=](const auto idx, const auto tidx) {
    to[idx][tidx] = from[idx][tidx];
  });
}

template<typename MDSfrom, typename MDSto>
void
  mds_add(const tyvi::mdgrid_work& w, MDSfrom&& from, MDSto&& to)
{

  if(from.extents() != to.extents()) {
    throw std::runtime_error { "Can not copy from/to different sized mdspan." };
  }

  w.for_each_index(from, [=](const auto idx, const auto tidx) {
    to[idx][tidx] = to[idx][tidx] + from[idx][tidx];
  });
}

}  // namespace

namespace emf {

YeeLattice::YeeLattice(const YeeLatticeCtorArgs args) :
  halo_size_ { args.halo_size },
  extents_wout_halo_ { args.Nx, args.Ny, args.Nz },
  E_(
    args.Nx + 2u * halo_size_,
    args.Ny + 2u * halo_size_,
    args.Nz + 2u * halo_size_),
  B_(
    args.Nx + 2u * halo_size_,
    args.Ny + 2u * halo_size_,
    args.Nz + 2u * halo_size_),
  J_(args.Nx + 2u * halo_size_, args.Ny + 2u * halo_size_, args.Nz + 2u * halo_size_)
{

  auto less_than_halo = [this](const auto x) { return x < halo_size_; };

  if(std::ranges::any_of(extents_wout_halo_, less_than_halo)) {
    auto msg = std::stringstream {};

    msg << "Yee Lattice extents (" << extents_wout_halo_[0] << ", "
        << extents_wout_halo_[1] << ", " << extents_wout_halo_[2]
        << ") are assumed to be at least halo size: " << halo_size_;
    throw std::runtime_error { msg.str() };
  }
}

std::array<std::size_t, 3>
  YeeLattice::extents_wout_halo() const
{
  return extents_wout_halo_;
}
std::array<std::size_t, 3>
  YeeLattice::extents_with_halo() const
{
  auto e        = extents_wout_halo();
  auto add_halo = [this](const auto x) { return x + 2uz * halo_size(); };
  std::ranges::transform(e, e.begin(), add_halo);
  return e;
}
std::size_t
  YeeLattice::halo_size() const
{
  return halo_size_;
}

YeeLattice::YeeLatticeHostCopy
  YeeLattice::get_EBJ()
{
  auto wE = tyvi::mdgrid_work {};
  auto wB = tyvi::mdgrid_work {};
  auto wJ = tyvi::mdgrid_work {};

  wE.sync_to_staging(E_);
  wB.sync_to_staging(B_);
  wJ.sync_to_staging(J_);

  const auto [Emds, Bmds, Jmds] = std::tuple { nonhalo_submds(E_.staging_mds()),
                                               nonhalo_submds(B_.staging_mds()),
                                               nonhalo_submds(J_.staging_mds()) };

  const auto [a, b, c] = extents_wout_halo();
  auto host_buffer     = YeeLatticeHostCopy(a, b, c);
  const auto mds       = host_buffer.mds();

  tyvi::when_all(wE, wB, wJ);
  wE.wait();

  for(const auto idx: tyvi::sstd::index_space(mds)) {
    mds[idx][] = YeeLatticeFieldsAtPoint { .Ex = Emds[idx][0],
                                           .Ey = Emds[idx][1],
                                           .Ez = Emds[idx][2],
                                           .Bx = Bmds[idx][0],
                                           .By = Bmds[idx][1],
                                           .Bz = Bmds[idx][2],
                                           .Jx = Jmds[idx][0],
                                           .Jy = Jmds[idx][1],
                                           .Jz = Jmds[idx][2] };
  }

  return host_buffer;
}

YeeLattice::YeeLatticeHostCopy
  YeeLattice::get_EBJ_with_halo()
{
  auto wE = tyvi::mdgrid_work {};
  auto wB = tyvi::mdgrid_work {};
  auto wJ = tyvi::mdgrid_work {};

  wE.sync_to_staging(E_);
  wB.sync_to_staging(B_);
  wJ.sync_to_staging(J_);

  const auto [Emds, Bmds, Jmds] =
    std::tuple { E_.staging_mds(), B_.staging_mds(), J_.staging_mds() };

  const auto [a, b, c] = extents_with_halo();
  auto host_buffer     = YeeLatticeHostCopy(a, b, c);
  const auto mds       = host_buffer.mds();

  tyvi::when_all(wE, wB, wJ);
  wE.wait();

  for(const auto idx: tyvi::sstd::index_space(mds)) {
    mds[idx][] = YeeLatticeFieldsAtPoint { .Ex = Emds[idx][0],
                                           .Ey = Emds[idx][1],
                                           .Ez = Emds[idx][2],
                                           .Bx = Bmds[idx][0],
                                           .By = Bmds[idx][1],
                                           .Bz = Bmds[idx][2],
                                           .Jx = Jmds[idx][0],
                                           .Jy = Jmds[idx][1],
                                           .Jz = Jmds[idx][2] };
  }

  return host_buffer;
}

void
  YeeLattice::add_current()
{
  const auto w = tyvi::mdgrid_work {};
  add_current(w);
  w.wait();
}

void
  YeeLattice::add_current(const tyvi::mdgrid_work& w)
{
  const auto Emds = nonhalo_submds(E_.mds());
  const auto Jmds = nonhalo_submds(J_.mds());

  w.for_each_index(Emds, [=](const auto idx, const auto tidx) {
    Emds[idx][tidx] = Emds[idx][tidx] - Jmds[idx][tidx];
  });
}

void
  YeeLattice::set_E_in_subregion(const dir_type dir, const YeeLattice& other)
{
  auto w = tyvi::mdgrid_work {};
  set_E_in_subregion(w, dir, other);
  w.wait();
}

void
  YeeLattice::set_B_in_subregion(const dir_type dir, const YeeLattice& other)
{
  auto w = tyvi::mdgrid_work {};
  set_B_in_subregion(w, dir, other);
  w.wait();
}

void
  YeeLattice::set_J_in_subregion(const dir_type dir, const YeeLattice& other)
{
  auto w = tyvi::mdgrid_work {};
  set_J_in_subregion(w, dir, other);
  w.wait();
}

void
  YeeLattice::set_E_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const YeeLattice& other)
{
  const auto my_Emds_region    = this->subregion(dir, this->E_.mds());
  const auto other_Emds_region = other.corresponding_subregion(dir, other.E_.mds());

  mds_copy(w, other_Emds_region, my_Emds_region);
}

void
  YeeLattice::set_B_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const YeeLattice& other)
{
  const auto my_Bmds_region    = this->subregion(dir, this->B_.mds());
  const auto other_Bmds_region = other.corresponding_subregion(dir, other.B_.mds());

  mds_copy(w, other_Bmds_region, my_Bmds_region);
}

void
  YeeLattice::set_J_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const YeeLattice& other)
{
  const auto my_Jmds_region    = this->subregion(dir, this->J_.mds());
  const auto other_Jmds_region = other.corresponding_subregion(dir, other.J_.mds());

  mds_copy(w, other_Jmds_region, my_Jmds_region);
}

void
  YeeLattice::add_to_J_from_subregion(const dir_type dir, const YeeLattice& other)
{
  auto w = tyvi::mdgrid_work {};
  add_to_J_from_subregion(w, dir, other);
  w.wait();
}

void
  YeeLattice::add_to_J_from_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const YeeLattice& other)
{
  const auto idir = invert_dir(dir);

  const auto my_Jmds_region    = this->corresponding_subregion(idir, this->J_.mds());
  const auto other_Jmds_region = other.subregion(idir, other.J_.mds());

  mds_add(w, other_Jmds_region, my_Jmds_region);
}

void
  YeeLattice::zero_transverse_E_behind_x(const std::size_t ix_wall)
{
  if(ix_wall == 0) return;

  const auto [nx, ny, nz] = extents_with_halo();
  const auto ix_end        = std::min(ix_wall + 1, nx);

  const auto Emds  = E_.mds();
  const auto region = std::submdspan(
    Emds,
    std::tuple { std::size_t { 0 }, ix_end },
    std::tuple { std::size_t { 0 }, ny },
    std::tuple { std::size_t { 0 }, nz });

  tyvi::mdgrid_work {}
    .for_each_index(
      region,
      [=](const auto idx) {
        region[idx][1] = 0;  // E_y = 0
        region[idx][2] = 0;  // E_z = 0
      })
    .wait();
}

void
  YeeLattice::clear_current()
{
  const tyvi::mdgrid_work w {};
  clear_current(w);
  w.wait();
}

void
  YeeLattice::clear_current(const tyvi::mdgrid_work& w)
{
  const auto Jmds = J_.mds();
  w.for_each_index(Jmds, [=](const auto idx, const auto tidx) { Jmds[idx][tidx] = 0; });
}

void
  YeeLattice::deposit_current(const CurrentContributions& contributions)
{
  const tyvi::mdgrid_work w {};
  deposit_current(w, contributions);
  w.wait();
}

void
  YeeLattice::deposit_current(
    const tyvi::mdgrid_work& w,
    const CurrentContributions& contributions)
{
  const auto b = thrust::make_zip_iterator(
    contributions.locations.begin(),
    contributions.currents.begin());
  const auto e = thrust::make_zip_iterator(
    contributions.locations.end(),
    contributions.currents.end());

  thrust::for_each(w.on_this(), b, e, [Jmds = J_.mds()](const auto& i) {
    const auto idx = thrust::get<0>(i);
    const auto J   = thrust::get<1>(i);
    Jmds[idx][0]   = Jmds[idx][0] + J[0];
    Jmds[idx][1]   = Jmds[idx][1] + J[1];
    Jmds[idx][2]   = Jmds[idx][2] + J[2];
  });
}

void
  YeeLattice::deposit_current(const runko::VecGrid<value_type>& depJ)
{
  const tyvi::mdgrid_work w {};
  deposit_current(w, depJ);
  w.wait();
}

void
  YeeLattice::deposit_current(
    const tyvi::mdgrid_work& w,
    const runko::VecGrid<value_type>& depJ)
{

  if(depJ.extents() != this->extents_with_halo()) {
    throw std::runtime_error {
      "emf::Tile::deposit_current: given current grid has incorrect extents!"
    };
  }
  const auto Jmds = this->J_.mds();
  w.for_each_index(Jmds, [=, depJmds = depJ.mds()](const auto idx, const auto tidx) {
    Jmds[idx][tidx] = Jmds[idx][tidx] + depJmds[idx][tidx];
  });
}

YeeLattice::VecGridMDS::mapping_type
  YeeLattice::grid_mapping_with_halo() const noexcept
{
  return this->E_.mds().mapping();
}

double
  YeeLattice::total_energy_B() const
{
  auto w = tyvi::mdgrid_work {};

  namespace rn           = std::ranges;
  const auto Bmds        = this->B_.mds();
  const auto index_space = tyvi::sstd::index_space(Bmds);

  const auto index_to_B2 = [=](const auto idx) {
    const auto Bx = Bmds[idx][0];
    const auto By = Bmds[idx][1];
    const auto Bz = Bmds[idx][2];
    return Bx * Bx + By * By + Bz * Bz;
  };
  const auto B2_iterator_begin =
    thrust::make_transform_iterator(index_space.begin(), index_to_B2);
  const auto B2_iterator_end = rn::next(B2_iterator_begin, rn::size(index_space));

  static constexpr auto norm = 8.0 * std::numbers::pi_v<double>;
  return thrust::reduce(w.on_this(), B2_iterator_begin, B2_iterator_end) / norm;
}


double
  YeeLattice::total_energy_E() const
{
  auto w = tyvi::mdgrid_work {};

  namespace rn           = std::ranges;
  const auto Emds        = this->E_.mds();
  const auto index_space = tyvi::sstd::index_space(Emds);

  const auto index_to_E2 = [=](const auto idx) {
    const auto Ex = Emds[idx][0];
    const auto Ey = Emds[idx][1];
    const auto Ez = Emds[idx][2];
    return Ex * Ex + Ey * Ey + Ez * Ez;
  };
  const auto E2_iterator_begin =
    thrust::make_transform_iterator(index_space.begin(), index_to_E2);
  const auto E2_iterator_end = rn::next(E2_iterator_begin, rn::size(index_space));

  static constexpr auto norm = 8.0 * std::numbers::pi_v<double>;
  return thrust::reduce(w.on_this(), E2_iterator_begin, E2_iterator_end) / norm;
}

}  // namespace emf
