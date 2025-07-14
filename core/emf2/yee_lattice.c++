#include "core/emf2/yee_lattice.h"

#include <algorithm>
#include <ranges>
#include <sstream>
#include <string>

namespace emf2 {

YeeLattice::YeeLattice(const YeeLatticeCtorArgs args) :
  halo_size_ { args.halo_size },
  extents_wout_halo_ { args.Nx, args.Ny, args.Nz },
  E_(
    args.Nx + 2uz * halo_size_,
    args.Ny + 2uz * halo_size_,
    args.Nz + 2uz * halo_size_),
  B_(
    args.Nx + 2uz * halo_size_,
    args.Ny + 2uz * halo_size_,
    args.Nz + 2uz * halo_size_),
  J_(args.Nx + 2uz * halo_size_, args.Ny + 2uz * halo_size_, args.Nz + 2uz * halo_size_)
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

  auto wE1 = wE.sync_to_staging(E_);
  auto wB1 = wB.sync_to_staging(B_);
  auto wJ1 = wJ.sync_to_staging(J_);

  const auto [Emds, Bmds, Jmds] = std::tuple { nonhalo_submds(E_.staging_mds()),
                                               nonhalo_submds(B_.staging_mds()),
                                               nonhalo_submds(J_.staging_mds()) };

  const auto [a, b, c] = extents_wout_halo();
  auto host_buffer     = YeeLatticeHostCopy(a, b, c);
  const auto mds       = host_buffer.mds();

  tyvi::when_all(wE1, wB1, wJ1).wait();

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

  auto wE1 = wE.sync_to_staging(E_);
  auto wB1 = wB.sync_to_staging(B_);
  auto wJ1 = wJ.sync_to_staging(J_);

  const auto [Emds, Bmds, Jmds] =
    std::tuple { E_.staging_mds(), B_.staging_mds(), J_.staging_mds() };

  const auto [a, b, c] = extents_with_halo();
  auto host_buffer     = YeeLatticeHostCopy(a, b, c);
  const auto mds       = host_buffer.mds();

  tyvi::when_all(wE1, wB1, wJ1).wait();

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
  YeeLattice::add_J_to_E()
{
  auto w0 = tyvi::mdgrid_work {};

  const auto Emds = nonhalo_submds(E_.mds());
  const auto Jmds = nonhalo_submds(J_.mds());

  auto w1 = w0.for_each_index(Emds, [=](const auto idx, const auto tidx) {
    Emds[idx][tidx] = Emds[idx][tidx] + Jmds[idx][tidx];
  });

  w1.wait();
}

}  // namespace emf2
