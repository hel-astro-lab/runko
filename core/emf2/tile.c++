#include "core/emf2/tile.h"

#include "tyvi/mdspan.h"

#include <iostream>

namespace emf2 {

template<std::size_t D>
Tile<D>::Tile(
  const std::array<std::size_t, 3> tile_indices,
  const toolbox::ConfigParser& p) :
  corgi::Tile<D>(),
  yee_lattice_extents_wout_halo_ { p.get_or_throw<std::size_t>("NxMesh"),
                                   p.get_or_throw<std::size_t>("NyMesh"),
                                   p.get_or_throw<std::size_t>("NzMesh") },
  yee_lattice_(
    yee_lattice_extents_wout_halo_[0] + 2 * halo_length,
    yee_lattice_extents_wout_halo_[1] + 2 * halo_length,
    yee_lattice_extents_wout_halo_[2] + 2 * halo_length),
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
  const auto Lx = yee_lattice_extents_wout_halo_[0];
  const auto Ly = yee_lattice_extents_wout_halo_[1];
  const auto Lz = yee_lattice_extents_wout_halo_[2];

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
    const auto x_coeff =
      static_cast<double>(i) / static_cast<double>(yee_lattice_extents_wout_halo_[0]);
    const auto y_coeff =
      static_cast<double>(j) / static_cast<double>(yee_lattice_extents_wout_halo_[1]);
    const auto z_coeff =
      static_cast<double>(k) / static_cast<double>(yee_lattice_extents_wout_halo_[2]);

    return std::tuple<double, double, double> {
      static_cast<double>(this->mins[0]) + x_coeff * Lx,
      static_cast<double>(this->mins[1]) + y_coeff * Ly,
      static_cast<double>(this->mins[2]) + z_coeff * Lz
    };
  };

  const auto [Emds, Bmds, Jmds] = yee_lattice_staging_mds_wout_halo();

  for(const auto idx: tyvi::sstd::index_space(Emds)) {
    const auto [x, y, z] = global_coordinates(idx[0], idx[1], idx[2]);
    const auto ex        = E(x + 0.5, y, z);
    const auto ey        = E(x, y + 0.5, z);
    const auto ez        = E(x, y, z + 0.5);
    const auto bx        = B(x, y + 0.5, z + 0.5);
    const auto by        = B(x + 0.5, y, z + 0.5);
    const auto bz        = B(x + 0.5, y + 0.5, z);
    const auto jx        = J(x + 0.5, y, z);
    const auto jy        = J(x, y + 0.5, z);
    const auto jz        = J(x, y, z + 0.5);

    Emds[idx][0] = std::get<0>(ex);
    Emds[idx][1] = std::get<1>(ey);
    Emds[idx][2] = std::get<2>(ez);

    Bmds[idx][0] = std::get<0>(bx);
    Bmds[idx][1] = std::get<1>(by);
    Bmds[idx][2] = std::get<2>(bz);

    Jmds[idx][0] = std::get<0>(jx);
    Jmds[idx][1] = std::get<1>(jy);
    Jmds[idx][2] = std::get<2>(jz);
  }

  auto wE = tyvi::mdgrid_work {};
  auto wB = tyvi::mdgrid_work {};
  auto wJ = tyvi::mdgrid_work {};

  auto wE1 = wE.sync_from_staging(yee_lattice_.E);
  auto wB1 = wB.sync_from_staging(yee_lattice_.B);
  auto wJ1 = wJ.sync_from_staging(yee_lattice_.J);

  tyvi::when_all(wE1, wB1, wJ1).wait();
}

template<std::size_t D>
Tile<D>::host_EBJ_grids
  Tile<D>::get_EBJ()
{

  auto wE = tyvi::mdgrid_work {};
  auto wB = tyvi::mdgrid_work {};
  auto wJ = tyvi::mdgrid_work {};

  auto wE1 = wE.sync_to_staging(yee_lattice_.E);
  auto wB1 = wB.sync_to_staging(yee_lattice_.B);
  auto wJ1 = wJ.sync_to_staging(yee_lattice_.J);

  const auto [Emds, Bmds, Jmds] = yee_lattice_staging_mds_wout_halo();

  auto host_buffer =
    host_EBJ_grids { .E = YeeLattice::host_vec_grid(yee_lattice_extents_wout_halo_),
                     .B = YeeLattice::host_vec_grid(yee_lattice_extents_wout_halo_),
                     .J = YeeLattice::host_vec_grid(yee_lattice_extents_wout_halo_) };

  const auto hostEmds = host_buffer.E.mds();
  const auto hostBmds = host_buffer.B.mds();
  const auto hostJmds = host_buffer.J.mds();

  tyvi::when_all(wE1, wB1, wJ1).wait();

  for(const auto idx: tyvi::sstd::index_space(hostEmds)) {
    for(const auto tidx: tyvi::sstd::index_space(hostEmds[idx])) {
      hostEmds[idx][tidx] = Emds[idx][tidx];
      hostBmds[idx][tidx] = Bmds[idx][tidx];
      hostJmds[idx][tidx] = Jmds[idx][tidx];
    }
  }

  return host_buffer;
}

template<std::size_t D>
std::array<std::size_t, 3>
  Tile<D>::extents_wout_halo() const
{
  return yee_lattice_extents_wout_halo_;
}

}  // namespace emf2

template class emf2::Tile<3>;
