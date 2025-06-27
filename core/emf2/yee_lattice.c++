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

}  // namespace emf2
