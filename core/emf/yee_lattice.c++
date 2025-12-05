#include "core/emf/yee_lattice.h"

#include "thrust/execution_policy.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/tuple.h"

#include <algorithm>
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
  extents_wout_halo_ { args.Nx, args.Ny, args.Nz },
  E_(args.Nx + 2uz * halo_size, args.Ny + 2uz * halo_size, args.Nz + 2uz * halo_size),
  B_(args.Nx + 2uz * halo_size, args.Ny + 2uz * halo_size, args.Nz + 2uz * halo_size),
  J_(args.Nx + 2uz * halo_size, args.Ny + 2uz * halo_size, args.Nz + 2uz * halo_size)
{

  auto less_than_halo = [this](const auto x) { return x < halo_size; };

  if(std::ranges::any_of(extents_wout_halo_, less_than_halo)) {
    auto msg = std::stringstream {};

    msg << "Yee Lattice extents (" << extents_wout_halo_[0] << ", "
        << extents_wout_halo_[1] << ", " << extents_wout_halo_[2]
        << ") are assumed to be at least halo size: " << halo_size;
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
  auto add_halo = [this](const auto x) { return x + 2uz * halo_size; };
  std::ranges::transform(e, e.begin(), add_halo);
  return e;
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
  YeeLattice::subtract_J_from_E()
{
  const auto w = tyvi::mdgrid_work {};
  subtract_J_from_E(w);
  w.wait();
}

void
  YeeLattice::subtract_J_from_E(const tyvi::mdgrid_work& w)
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


void
  YeeLattice::set_E_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const toolbox::hollow_grid<value_type, 3, halo_size>& other)
{
  const auto my_Emds_region = this->subregion(dir, this->E_.mds());
  const auto acc            = other.get_accessor();
  const auto base =
    std::array { dir[0] == -1 ? extents_wout_halo_[0] - halo_size : 0uz,
                 dir[1] == -1 ? extents_wout_halo_[1] - halo_size : 0uz,
                 dir[2] == -1 ? extents_wout_halo_[2] - halo_size : 0uz };

  w.for_each_index(my_Emds_region, [=](const auto idx, const auto tidx) {
    const auto i              = std::array<std::size_t, 3> { idx[0] + base[0],
                                                             idx[1] + base[1],
                                                             idx[2] + base[2] };
    const auto n              = acc.offset(i, tidx[0]);
    my_Emds_region[idx][tidx] = acc.data[n];
  });
}

void
  YeeLattice::set_B_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const toolbox::hollow_grid<value_type, 3, halo_size>& other)
{
  const auto my_Bmds_region = this->subregion(dir, this->B_.mds());
  const auto acc            = other.get_accessor();
  const auto base =
    std::array { dir[0] == -1 ? extents_wout_halo_[0] - halo_size : 0uz,
                 dir[1] == -1 ? extents_wout_halo_[1] - halo_size : 0uz,
                 dir[2] == -1 ? extents_wout_halo_[2] - halo_size : 0uz };

  w.for_each_index(my_Bmds_region, [=](const auto idx, const auto tidx) {
    const auto i              = std::array<std::size_t, 3> { idx[0] + base[0],
                                                             idx[1] + base[1],
                                                             idx[2] + base[2] };
    const auto n              = acc.offset(i, tidx[0]);
    my_Bmds_region[idx][tidx] = acc.data[n];
  });
}

void
  YeeLattice::set_J_in_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const toolbox::hollow_grid<value_type, 3, 2 * halo_size>& other)
{
  const auto my_Jmds_region = this->subregion(dir, this->J_.mds());
  const auto acc            = other.get_accessor();

  auto base = std::array { dir[0] == -1 ? extents_wout_halo_[0] : halo_size,
                           dir[1] == -1 ? extents_wout_halo_[1] : halo_size,
                           dir[2] == -1 ? extents_wout_halo_[2] : halo_size };

  w.for_each_index(my_Jmds_region, [=](const auto idx, const auto tidx) {
    const auto i              = std::array<std::size_t, 3> { idx[0] + base[0],
                                                             idx[1] + base[1],
                                                             idx[2] + base[2] };
    const auto n              = acc.offset(i, tidx[0]);
    my_Jmds_region[idx][tidx] = acc.data[n];
  });
}

void
  YeeLattice::add_to_J_from_subregion(
    const tyvi::mdgrid_work& w,
    const dir_type dir,
    const toolbox::hollow_grid<value_type, 3, 2 * halo_size>& other)
{
  const auto idir           = invert_dir(dir);
  const auto my_Jmds_region = this->corresponding_subregion(idir, this->J_.mds());
  const auto acc            = other.get_accessor();

  auto base = std::array<std::size_t, 3> {};
  for(const auto i: std::array { 0uz, 1uz, 2uz }) {
    if(dir[i] == 1) {
      base[i] = 0uz;
    } else if(dir[i] == 0) {
      base[i] = halo_size;
    } else {
      base[i] = halo_size + extents_wout_halo_[i];
    }
  }

  w.for_each_index(my_Jmds_region, [=](const auto idx, const auto tidx) {
    const auto i              = std::array<std::size_t, 3> { idx[0] + base[0],
                                                             idx[1] + base[1],
                                                             idx[2] + base[2] };
    const auto n              = acc.offset(i, tidx[0]);
    my_Jmds_region[idx][tidx] = my_Jmds_region[idx][tidx] + acc.data[n];
  });
}

}  // namespace emf
