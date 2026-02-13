#include "core/emf/tile.h"

#include "core/communication_common.h"
#include "core/mdgrid_common.h"
#include "tools/system.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <cmath>
#include <format>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace {
emf::FieldPropagator
  parse_field_propagator(const std::string_view p)
{
  if(p == "FDTD2") {
    return emf::FieldPropagator::FDTD2;
  } else {
    const auto msg = std::format("{} is not supported field propagator.", p);
    throw std::runtime_error { msg };
  }
}

emf::CurrentFilter
  parse_current_filter(const std::string_view p)
{
  if(p == "binomial2") {
    return emf::CurrentFilter::binomial2;
  } else {
    const auto msg = std::format("{} is not supported current filter.", p);
    throw std::runtime_error { msg };
  }
}
}  // namespace

namespace emf {

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
  if(const auto cfilter = p.get<std::string>("current_filter")) {
    current_filter_ = parse_current_filter(cfilter.value());
  }


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
      "failed emf::tile::set_fields precondition: corgi::tile::{mins,maxs} aren't set"
    };
  }

  const auto global_coordinates = global_coordinate_map();

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
void
  Tile<D>::batch_set_EBJ(
    batch_vector_field_function Ex,
    batch_vector_field_function Ey,
    batch_vector_field_function Ez,
    batch_vector_field_function Bx,
    batch_vector_field_function By,
    batch_vector_field_function Bz,
    batch_vector_field_function Jx,
    batch_vector_field_function Jy,
    batch_vector_field_function Jz)
{
  if(this->mins == this->maxs) {
    throw std::logic_error {
      "failed emf::tile::set_fields precondition: corgi::tile::{mins,maxs} aren't set"
    };
  }

  const auto global_coordinates = global_coordinate_map();

  const auto e = yee_lattice_.extents_wout_halo();

  auto x   = pybind11::array_t<double>(e);
  auto y   = pybind11::array_t<double>(e);
  auto z   = pybind11::array_t<double>(e);
  auto xp5 = pybind11::array_t<double>(e);
  auto yp5 = pybind11::array_t<double>(e);
  auto zp5 = pybind11::array_t<double>(e);

  auto xv   = x.template mutable_unchecked<3>();
  auto yv   = y.template mutable_unchecked<3>();
  auto zv   = z.template mutable_unchecked<3>();
  auto xp5v = xp5.template mutable_unchecked<3>();
  auto yp5v = yp5.template mutable_unchecked<3>();
  auto zp5v = zp5.template mutable_unchecked<3>();


  for(const auto [i, j, k]:
      tyvi::sstd::index_space(std::mdspan((int*)nullptr, e[0], e[1], e[2]))) {

    const auto [gx, gy, gz]       = global_coordinates(i, j, k);
    const auto [gxp5, gyp5, gzp5] = global_coordinates(i + 0.5, j + 0.5, k + 0.5);
    xv(i, j, k)                   = gx;
    yv(i, j, k)                   = gy;
    zv(i, j, k)                   = gz;
    xp5v(i, j, k)                 = gxp5;
    yp5v(i, j, k)                 = gyp5;
    zp5v(i, j, k)                 = gzp5;
  }

  const auto ex = Ex(xp5, y, z);
  const auto ey = Ey(x, yp5, z);
  const auto ez = Ez(x, y, zp5);

  const auto bx = Bx(x, yp5, zp5);
  const auto by = By(xp5, y, zp5);
  const auto bz = Bz(xp5, yp5, z);

  const auto jx = Jx(xp5, y, z);
  const auto jy = Jy(x, yp5, z);
  const auto jz = Jz(x, y, zp5);

  const auto assert_shape = [&](const auto& A) {
    if(A.shape(0) != e[0] or A.shape(1) != e[1] or A.shape(2) != e[2]) {
      throw std::runtime_error {
        "Batch field setter returned array with incorrect shape!"
      };
    }
  };

  assert_shape(ex);
  assert_shape(ey);
  assert_shape(ez);
  assert_shape(bx);
  assert_shape(by);
  assert_shape(bz);
  assert_shape(jx);
  assert_shape(jy);
  assert_shape(jz);

  const auto exv = ex.template unchecked<3>();
  const auto eyv = ey.template unchecked<3>();
  const auto ezv = ez.template unchecked<3>();
  const auto bxv = bx.template unchecked<3>();
  const auto byv = by.template unchecked<3>();
  const auto bzv = bz.template unchecked<3>();
  const auto jxv = jx.template unchecked<3>();
  const auto jyv = jy.template unchecked<3>();
  const auto jzv = jz.template unchecked<3>();


  auto f = [&](const std::size_t i, const std::size_t j, const std::size_t k) {
    return YeeLatticeFieldsAtPoint { .Ex = exv(i, j, k),
                                     .Ey = eyv(i, j, k),
                                     .Ez = ezv(i, j, k),
                                     .Bx = bxv(i, j, k),
                                     .By = byv(i, j, k),
                                     .Bz = bzv(i, j, k),
                                     .Jx = jxv(i, j, k),
                                     .Jy = jyv(i, j, k),
                                     .Jz = jzv(i, j, k) };
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
        "emf::Tile::push_half_b internal error: field_propagator_ not set."
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
        "emf::Tile::push_e internal error: field_propagator_ not set."
      };
  }
}

template<std::size_t D>
void
  Tile<D>::subtract_J_from_E()
{
  yee_lattice_.subtract_J_from_E();
}

template<std::size_t D>
void
  Tile<D>::filter_current()
{
  if(not current_filter_) {
    throw std::logic_error {
      "Trying to filter current without specifying `current_filter`!"
    };
  }

  const auto cf = current_filter_.value();
  switch(cf) {
    case emf::CurrentFilter::binomial2:
      this->yee_lattice_.filter_current_binomial2();
      return;
    default:
      throw std::logic_error {
        "emf::Tile::filter_current internal error: unregonized current filter."
      };
  }
}

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

  auto make_isend = [&](const auto s) {
    return comm.isend(dest, tag, s.data(), s.size());
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E: return { make_isend(yee_lattice_.span_E()) };
    case comm_mode::emf_B: return { make_isend(yee_lattice_.span_B()) };
    case comm_mode::emf_J: return { make_isend(yee_lattice_.span_J()) };
    default:
      throw std::logic_error { std::format(
        "emf::Tile::send_data does not support given communication mode: {}",
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

  auto make_irecv = [&](const auto s) {
    return comm.irecv(orig, tag, s.data(), s.size());
  };

  switch(static_cast<comm_mode>(mode)) {
    case comm_mode::emf_E: return { make_irecv(yee_lattice_.span_E()) };
    case comm_mode::emf_B: return { make_irecv(yee_lattice_.span_B()) };
    case comm_mode::emf_J: return { make_irecv(yee_lattice_.span_J()) };
    default:
      throw std::logic_error { std::format(
        "emf::Tile::recv_data does not support given communication mode: {}",
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
        "emf::Tile::pairwise_moore_communication assumes that the other tile is "
        "emf::Tile or its descendant. Orginal exception: {}",
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
    case comm_mode::emf_J:
      yee_lattice_.set_J_in_subregion(dir_to_other, other.yee_lattice_);
      break;
    case comm_mode::emf_J_exchange:
      yee_lattice_.add_to_J_from_subregion(dir_to_other, other.yee_lattice_);
      break;
    default:
      throw std::logic_error { std::format(
        "emf::Tile::pairwise_moore_communication does not support given communication "
        "mode: {}",
        mode) };
  }
};

template<std::size_t D>
void
  Tile<D>::register_antenna(emf::antenna_mode mode)
{
  this->antenna_modes_.push_back(std::move(mode));
};

template<std::size_t D>
void
  Tile<D>::deposit_antenna_current()
{
  using vec_list = runko::VecList<emf::YeeLattice::value_type>;

  const auto num_of_modes = this->antenna_modes_.size();
  auto A                  = vec_list(num_of_modes);
  auto K                  = vec_list(num_of_modes);

  const auto sA_mds = A.staging_mds();
  const auto sk_mds = K.staging_mds();

  for(const auto n: std::views::iota(0uz, num_of_modes)) {
    for(const auto i: std::views::iota(0uz, 3uz)) {
      sA_mds[n][i] =
        static_cast<emf::YeeLattice::value_type>(this->antenna_modes_[n].A[i]);
      sk_mds[n][i] =
        static_cast<emf::YeeLattice::value_type>(this->antenna_modes_[n].k[i]);
    }
  }

  const tyvi::mdgrid_work w {};

  w.sync_from_staging(A).sync_from_staging(K);

  auto vec_pot =
    runko::VecGrid<emf::YeeLattice::value_type>(this->yee_lattice_.extents_with_halo());

  const auto A_mds       = A.mds();
  const auto K_mds       = K.mds();
  const auto vec_pot_mds = vec_pot.mds();

  const auto gcmap = this->global_coordinate_map();

  w.for_each_index(vec_pot_mds, [=](const auto idx) {
    // For each is over all indices (including halo region)
    // but the global coordinate map is defined s.t. (0, 0, 0) is located at corner of
    // non-halo region.
    const auto i = static_cast<double>(idx[0]) - emf::Tile<D>::halo_size;
    const auto j = static_cast<double>(idx[1]) - emf::Tile<D>::halo_size;
    const auto k = static_cast<double>(idx[2]) - emf::Tile<D>::halo_size;

    using vec        = toolbox::Vec3<double>;
    const auto x_loc = vec(gcmap(i + 0.5, j, k));
    const auto y_loc = vec(gcmap(i, j + 0.5, k));
    const auto z_loc = vec(gcmap(i, j, k + 0.5));

    for(auto n = 0uz; n < num_of_modes; ++n) {
      vec_pot_mds[idx][0] = vec_pot_mds[idx][0] +
                            A_mds[n][0] * std::cos(toolbox::dot(x_loc, vec(K_mds[n])));
      vec_pot_mds[idx][1] = vec_pot_mds[idx][1] +
                            A_mds[n][1] * std::cos(toolbox::dot(y_loc, vec(K_mds[n])));
      vec_pot_mds[idx][2] = vec_pot_mds[idx][2] +
                            A_mds[n][2] * std::cos(toolbox::dot(z_loc, vec(K_mds[n])));
    }
  });

  // We have to calculate B = curl(vec_pot) only in non-halo region + one deep shell in
  // halo region.
  auto generated_B =
    runko::VecGrid<emf::YeeLattice::value_type>(this->yee_lattice_.extents_with_halo());

  const auto B_mds = generated_B.mds();

  const auto h   = this->halo_size;
  const auto hp1 = h + 1uz;

  const auto [ex, ey, ez] = this->extents_wout_halo();
  const auto i1           = std::tuple { h - 1uz, h + ex + 1uz };
  const auto j1           = std::tuple { h - 1uz, h + ey + 1uz };
  const auto k1           = std::tuple { h - 1uz, h + ez + 1uz };

  const auto i1p1 = std::tuple { h, h + ex + 2uz };
  const auto j1p1 = std::tuple { h, h + ey + 2uz };
  const auto k1p1 = std::tuple { h, h + ez + 2uz };

  // FIXME: unify this with emf FDTD2.
  auto curl = [&](
                const auto coeff,
                const auto out,
                const auto X,
                const auto Xip1,
                const auto Xjp1,
                const auto Xkp1) {
    w.for_each_index(out, [=](const auto idx) {
      const auto Dk = Xkp1[idx][1] - X[idx][1];
      const auto Dj = Xjp1[idx][2] - X[idx][2];
      out[idx][0]   = coeff * (Dj - Dk);
    });
    w.for_each_index(out, [=](const auto idx) {
      const auto Di = Xip1[idx][2] - X[idx][2];
      const auto Dk = Xkp1[idx][0] - X[idx][0];
      out[idx][1]   = coeff * (Dk - Di);
    });
    w.for_each_index(out, [=](const auto idx) {
      const auto Dj = Xjp1[idx][0] - X[idx][0];
      const auto Di = Xip1[idx][1] - X[idx][1];
      out[idx][2]   = coeff * (Di - Dj);
    });
  };

  curl(
    1,
    std::submdspan(B_mds, i1, j1, k1),
    std::submdspan(vec_pot_mds, i1, j1, k1),
    std::submdspan(vec_pot_mds, i1p1, j1, k1),
    std::submdspan(vec_pot_mds, i1, j1p1, k1),
    std::submdspan(vec_pot_mds, i1, j1, k1p1));


  // Now curl(B) in non-halo region.
  // We can reuse vec_pot container.

  const auto i   = std::tuple { h, h + ex };
  const auto j   = std::tuple { h, h + ey };
  const auto k   = std::tuple { h, h + ez };
  const auto im1 = std::tuple { h - 1uz, h + ex - 1uz };
  const auto jm1 = std::tuple { h - 1uz, h + ey - 1uz };
  const auto km1 = std::tuple { h - 1uz, h + ez - 1uz };

  // Here we negate cfl, as we put (im1, jm1, km1) instead of (ip1, jp1, kp1).
  curl(
    -this->cfl_,
    std::submdspan(vec_pot_mds, i, j, k),
    std::submdspan(B_mds, i, j, k),
    std::submdspan(B_mds, im1, j, k),
    std::submdspan(B_mds, i, jm1, k),
    std::submdspan(B_mds, i, j, km1));

  this->yee_lattice_.deposit_current(w, vec_pot);
  w.wait();
};


}  // namespace emf

template class emf::Tile<3>;
