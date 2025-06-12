#include "core/emf2/tile.h"

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
    yee_lattice_extents_[0],
    yee_lattice_extents_[1],
    yee_lattice_extents_[2]),
  cfl_ { p.get<double>("cfl").value() }
{
}

}  // namespace emf2

template class emf2::Tile<3>;
