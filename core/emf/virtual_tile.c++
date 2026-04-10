#include "core/emf/virtual_tile.h"

#include "core/communication_common.h"
#include "thrust/host_vector.h"
#include "tools/system.h"

namespace emf {

template<std::size_t D>
VirtualTile<D>::VirtualTile(
  const std::array<std::size_t, 3> tile_indices,
  const toolbox::ConfigParser& p) :
  corgi::Tile<D>(),
  extents_wout_halo_ { p.get_or_throw<std::size_t>("NxMesh"),
                       p.get_or_throw<std::size_t>("NyMesh"),
                       p.get_or_throw<std::size_t>("NzMesh") }
{

  E_ = hollow_grid_EB(extents_wout_halo());
  B_ = hollow_grid_EB(extents_wout_halo());
  J_ = hollow_grid_J(extents_with_halo());

  auto one_or_throw = [](const double d) {
    if(d != 1.0) {
      throw std::logic_error {
        "emf tile does not support d{x,y,z} values other than 1."
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

  const auto Nx = p.get_or_throw<std::size_t>("Nx");
  const auto Ny = p.get_or_throw<std::size_t>("Ny");
  const auto Nz = p.get_or_throw<std::size_t>("Nz");

  if(Nx <= i or Ny <= j or Nz <= k) {
    throw std::runtime_error { "Trying to create tile outside of configured grid." };
  }

  this->index = { tile_indices[0], tile_indices[1], tile_indices[2] };

  // Tile size.
  const auto [Lx, Ly, Lz] = extents_wout_halo_;

  this->set_tile_mins({ xmin + i * Lx, ymin + j * Ly, zmin + k * Lz });
  this->set_tile_maxs(
    { xmin + (i + 1uz) * Lx, ymin + (j + 1uz) * Ly, zmin + (k + 1uz) * Lz });
}

template<std::size_t D>
std::array<std::size_t, 3>
  VirtualTile<D>::extents_wout_halo() const
{
  return this->extents_wout_halo_;
}

template<std::size_t D>
std::array<std::size_t, 3>
  VirtualTile<D>::extents_with_halo() const
{

  return { this->extents_wout_halo_[0] + 2 * halo_size,
           this->extents_wout_halo_[1] + 2 * halo_size,
           this->extents_wout_halo_[2] + 2 * halo_size };
}

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  VirtualTile<D>::recv_data(
    mpi4cpp::mpi::communicator& comm,
    const int orig,
    const int mode,
    const int tag)
{
  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error {
      "Non gpu aware MPI communication is not yet implemented."
    };
  }

  using runko::comm_mode;

  auto make_irecv = [&](const auto s) {
    return comm.irecv(orig, tag, s.data(), s.size());
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E: return { make_irecv(E_.span()) };
    case comm_mode::emf_B: return { make_irecv(B_.span()) };
    case comm_mode::emf_J: return { make_irecv(J_.span()) };
    default:
      throw std::logic_error { std::format(
        "emf::Tile::recv_data does not support given communication mode: {}",
        mode) };
  }
}
}  // namespace emf

template class emf::VirtualTile<3>;
