#include "core/emf2/tile.h"

#include "tyvi/mdspan.h"

#include <iostream>

namespace emf2 {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_indices,
  const toolbox::ConfigParser& p) :
  corgi::Tile<D>(),
  yee_lattice_(
    YeeLatticeCtorArgs { .halo_size = halo_size,
                         .Nx        = p.get_or_throw<std::size_t>("NxMesh"),
                         .Ny        = p.get_or_throw<std::size_t>("NyMesh"),
                         .Nz        = p.get_or_throw<std::size_t>("NzMesh")

    }),
  cfl_ { p.get_or_throw<double>("cfl") }
{
  const auto Nx = p.get_or_throw<std::size_t>("Nx");
  const auto Ny = p.get_or_throw<std::size_t>("Ny");
  const auto Nz = p.get_or_throw<std::size_t>("Nz");

  auto one_or_throw = [](const double d) {
    if(d != 1.0) {
      throw std::logic_error {
        "emf2 tile does not support d{x,y,z} values other than 1."
      };
    }
    return std::optional { d };
  };

  std::ignore = p.get<double>("dx").and_then(one_or_throw);
  std::ignore = p.get<double>("dy").and_then(one_or_throw);
  std::ignore = p.get<double>("dz").and_then(one_or_throw);

  const auto xmin = p.get_or_throw<double>("xmin");
  const auto ymin = p.get_or_throw<double>("ymin");
  const auto zmin = p.get_or_throw<double>("zmin");

  const auto [i, j, k] = tile_indices;

  // Tile size.
  const auto [Lx, Ly, Lz] = yee_lattice_.extents_wout_halo();

  this->set_tile_mins({ xmin + i * Lx, ymin + j * Ly, zmin + k * Lz });
  this->set_tile_maxs(
    { xmin + (i + 1uz) * Lx, ymin + (j + 1uz) * Ly, zmin + (k + 1uz) * Lz });
}

template<std::size_t D>
void
  Tile<D>::set_EBJ(
    vector_field_function E,
    vector_field_function B,
    vector_field_function J)
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
    const auto e       = yee_lattice_.extents_wout_halo();
    const auto x_coeff = static_cast<double>(i) / static_cast<double>(e[0]);
    const auto y_coeff = static_cast<double>(j) / static_cast<double>(e[1]);
    const auto z_coeff = static_cast<double>(k) / static_cast<double>(e[2]);

    return std::tuple<double, double, double> {
      static_cast<double>(this->mins[0]) + x_coeff * Lx,
      static_cast<double>(this->mins[1]) + y_coeff * Ly,
      static_cast<double>(this->mins[2]) + z_coeff * Lz
    };
  };

  const auto [Emds, Bmds, Jmds] = yee_lattice_.staging_mds_wout_halo();

  auto f = [&](const std::size_t i, const std::size_t j, const std::size_t k) {
    const auto [x, y, z] = global_coordinates(i, j, k);
    const auto ex        = E(x + 0.5, y, z);
    const auto ey        = E(x, y + 0.5, z);
    const auto ez        = E(x, y, z + 0.5);
    const auto bx        = B(x, y + 0.5, z + 0.5);
    const auto by        = B(x + 0.5, y, z + 0.5);
    const auto bz        = B(x + 0.5, y + 0.5, z);
    const auto jx        = J(x + 0.5, y, z);
    const auto jy        = J(x, y + 0.5, z);
    const auto jz        = J(x, y, z + 0.5);

    return YeeLatticeFieldsAtPoint { .Ex = std::get<0>(ex),
                                     .Ey = std::get<1>(ey),
                                     .Ez = std::get<2>(ez),
                                     .Bx = std::get<0>(bx),
                                     .By = std::get<1>(by),
                                     .Bz = std::get<2>(bz),
                                     .Jx = std::get<0>(jx),
                                     .Jy = std::get<1>(jy),
                                     .Jz = std::get<2>(jz) };
  };

  yee_lattice_.set_EBJ(f);
}

template<std::size_t D>
YeeLattice::YeeLatticeHostCopy
  Tile<D>::get_EBJ()
{
  return yee_lattice_.get_EBJ();
}

template<std::size_t D>
std::array<std::size_t, 3>
  Tile<D>::extents_wout_halo() const
{
  return yee_lattice_.extents_wout_halo();
}

}  // namespace emf2

template class emf2::Tile<3>;
