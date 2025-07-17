#include "core/emf2/tile.h"

#include "core/communication_common.h"
#include "tools/system.h"
#include "tyvi/mdspan.h"

#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace {
emf2::FieldPropagator
  parse_field_propagator(const std::string_view p)
{
  if(p == "FDTD2") {
    return emf2::FieldPropagator::FDTD2;
  } else {
    const auto msg = std::format("{} is not supported field propagator.", p);
    throw std::runtime_error { msg };
  }
}
}  // namespace

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
  cfl_ { p.get_or_throw<double>("cfl") },
  field_propagator_ { parse_field_propagator(
    p.get_or_throw<std::string>("field_propagator")) }
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

  /// FIXME: unify global coordinates from here and pic2::Tile.
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
YeeLattice::YeeLatticeHostCopy
  Tile<D>::get_EBJ_with_halo()
{
  return yee_lattice_.get_EBJ_with_halo();
}

template<std::size_t D>
std::array<std::size_t, 3>
  Tile<D>::extents_wout_halo() const
{
  return yee_lattice_.extents_wout_halo();
}

template<std::size_t D>
void
  Tile<D>::push_half_b()
{
  switch(field_propagator_) {
    case FieldPropagator::FDTD2: yee_lattice_.push_b_FDTD2(cfl_ / 2); break;
    default:
      throw std::logic_error {
        "emf2::Tile::push_half_b internal error: field_propagator_ not set."
      };
  }
}

template<std::size_t D>
void
  Tile<D>::push_e()
{
  switch(field_propagator_) {
    case FieldPropagator::FDTD2: yee_lattice_.push_e_FDTD2(cfl_); break;
    default:
      throw std::logic_error {
        "emf2::Tile::push_e internal error: field_propagator_ not set."
      };
  }
}

template<std::size_t D>
void
  Tile<D>::deposit_current()
{
  yee_lattice_.add_J_to_E();
}

namespace {
template<runko::comm_mode m>
constexpr int
  specialized_tag(const int tag)
{

  // cray-mpich supports maximum of 2^22-1 tag value.
  if(tag / 3 + 1 >= std::pow(2, 22)) {
    throw std::runtime_error { "emf2::Tile does not support tags >= 2^22 / 3." };
  }

  static constexpr auto comm_ordinal = [] {
    using runko::comm_mode;
    static_assert(
      m == comm_mode::emf_E or m == comm_mode::emf_B or m == comm_mode::emf_J);

    switch(m) {
      case comm_mode::emf_E: return 0;
      case comm_mode::emf_B: return 1;
      case comm_mode::emf_J: return 2;
    }
  }();

  return 3 * tag + comm_ordinal;
}
}  // namespace

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  Tile<D>::send_data(
    mpi4cpp::mpi::communicator& comm,
    const int dest,
    const int mode,
    const int tag)
{

  if(not toolbox::system_supports_gpu_aware_mpi()) {
    throw std::runtime_error {
      "Non gpu aware MPI communication is not yet implemented."
    };
  }

  using runko::comm_mode;

  auto make_isend = [&](const auto s, const int t) {
    return comm.isend(dest, t, s.data(), s.size());
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E:
      return {
        make_isend(yee_lattice_.span_E(), specialized_tag<comm_mode::emf_E>(tag))
      };
    case comm_mode::emf_B:
      return {
        make_isend(yee_lattice_.span_B(), specialized_tag<comm_mode::emf_B>(tag))
      };
    case comm_mode::emf_J:
      return {
        make_isend(yee_lattice_.span_J(), specialized_tag<comm_mode::emf_J>(tag))
      };
    default:
      throw std::logic_error { std::format(
        "emf2::Tile::send_data does not support given communication mode: {}",
        mode) };
  }
}

template<std::size_t D>
std::vector<mpi4cpp::mpi::request>
  Tile<D>::recv_data(
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

  auto make_irecv = [&](const auto s, const int t) {
    return comm.irecv(orig, t, s.data(), s.size());
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E:
      return {
        make_irecv(yee_lattice_.span_E(), specialized_tag<comm_mode::emf_E>(tag))
      };
    case comm_mode::emf_B:
      return {
        make_irecv(yee_lattice_.span_B(), specialized_tag<comm_mode::emf_B>(tag))
      };
    case comm_mode::emf_J:
      return {
        make_irecv(yee_lattice_.span_J(), specialized_tag<comm_mode::emf_J>(tag))
      };
    default:
      throw std::logic_error { std::format(
        "emf2::Tile::recv_data does not support given communication mode: {}",
        mode) };
  }
}

template<std::size_t D>
void
  Tile<D>::pairwise_moore_communication(
    const corgi::Tile<D>& other_base,
    const std::array<int, D> dir_to_other,
    const int mode)
{
  const Tile<D>& other = [&]() -> const Tile<D>& {
    try {
      return dynamic_cast<const Tile<D>&>(other_base);
    } catch(const std::bad_cast& ex) {
      throw std::runtime_error { std::format(
        "emf2::Tile::pairwise_moore_communication assumes that the other tile is "
        "emf2::Tile or its descendant. Orginal exception: {}",
        ex.what()) };
    }
  }();

  using runko::comm_mode;

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E:
      yee_lattice_.set_E_in_subregion(dir_to_other, other.yee_lattice_);
      break;
    case comm_mode::emf_B:
      yee_lattice_.set_B_in_subregion(dir_to_other, other.yee_lattice_);
      break;
    default:
      throw std::logic_error { std::format(
        "emf2::Tile::pairwise_moore_communication does not support given communication "
        "mode: {}",
        mode) };
  }
};

}  // namespace emf2

template class emf2::Tile<3>;
