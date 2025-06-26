#include "core/emf2/yee_lattice.h"

namespace emf2 {

YeeLattice::YeeLattice(
  const std::size_t Nx,
  const std::size_t Ny,
  const std::size_t Nz) :
  E(Nx, Ny, Nz),
  B(Nx, Ny, Nz),
  J(Nx, Ny, Nz)
{
}

}  // namespace emf2
