#include "core/emf2/tile.h"

#include "tyvi/mdspan.h"

#include <iostream>

namespace emf2 {

YeeLattice::YeeLattice(
  const std::size_t Nx,
  const std::size_t Ny,
  const std::size_t Nz) :
  E(Nx, Ny, Nz),
  B(Nx, Ny, Nz),
  rho(Nx, Ny, Nz),
  J(Nx, Ny, Nz)
{
}

template<std::size_t D>
Tile<D>::Tile(const toolbox::ConfigParser& p) :
  corgi::Tile<D>(),
  yee_lattice_extents_ { p.get<std::size_t>("Nx").value(),
                         p.get<std::size_t>("Ny").value(),
                         p.get<std::size_t>("Nz").value() },
  yee_lattice_(
    yee_lattice_extents_[0] + 2 * halo_length,
    yee_lattice_extents_[1] + 2 * halo_length,
    yee_lattice_extents_[2] + 2 * halo_length),
  cfl_ { p.get<double>("cfl").value() }
{
}

template<std::size_t D>
void
  Tile<D>::set_fields(vector_field E, vector_field B, scalar_field rho, vector_field J)
{
  if(this->mins == this->maxs) {
    throw std::logic_error {
      "failed emf2::tile::set_fields precondition: corgi::tile::{mins,maxs} aren't set"
    };
  }

  const auto Lx = static_cast<double>(this->maxs[0] - this->mins[0]);
  const auto Ly = static_cast<double>(this->maxs[1] - this->mins[1]);
  const auto Lz = static_cast<double>(this->maxs[2] - this->mins[2]);

  auto global_coordinates = [&](const auto i, const auto j, const auto k) {
    const auto x_coeff =
      static_cast<double>(i) / static_cast<double>(yee_lattice_extents_[0]);
    const auto y_coeff =
      static_cast<double>(i) / static_cast<double>(yee_lattice_extents_[1]);
    const auto z_coeff =
      static_cast<double>(i) / static_cast<double>(yee_lattice_extents_[2]);

    return std::tuple<double, double, double> {
      static_cast<double>(this->mins[0]) + x_coeff * Lx,
      static_cast<double>(this->mins[1]) + y_coeff * Ly,
      static_cast<double>(this->mins[2]) + z_coeff * Lz
    };
  };

  const auto [Emds, Bmds, rhomds, Jmds] = yee_lattice_staging_mds_no_halo();

  for(const auto idx: tyvi::sstd::index_space(Emds)) {
    const auto [x, y, z] = global_coordinates(idx[0], idx[1], idx[2]);
    const auto e         = E(x, y, z);
    const auto b         = B(x, y, z);
    const auto j         = J(x, y, z);

    rhomds[idx][] = rho(x, y, z);

    for(const auto tidx: tyvi::sstd::index_space(Emds[idx])) {
      Emds[idx][tidx] = e[tidx[0]];
      Bmds[idx][tidx] = b[tidx[0]];
      Jmds[idx][tidx] = j[tidx[0]];
    }
  }

  auto wE   = tyvi::mdgrid_work {};
  auto wB   = tyvi::mdgrid_work {};
  auto wRho = tyvi::mdgrid_work {};
  auto wJ   = tyvi::mdgrid_work {};

  auto wE1   = wE.sync_from_staging(yee_lattice_.E);
  auto wB1   = wB.sync_from_staging(yee_lattice_.B);
  auto wRho1 = wRho.sync_from_staging(yee_lattice_.rho);
  auto wJ1   = wJ.sync_from_staging(yee_lattice_.J);

  tyvi::when_all(wE1, wB1, wRho1, wJ1).wait();
}

}  // namespace emf2

template class emf2::Tile<3>;
