#include "core/emf/antenna.h"
#include "core/emf/tile.h"
#include "io/emf_average_field_energy_density.h"
#include "io/snapshots/hdf5_fields.h"
#include "io/snapshots/mpiio_fields.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

namespace {

namespace py = pybind11;

/// MVP implementation.
auto
  to_ndarrays(emf::YeeLattice::YeeLatticeHostCopy& lattice)
{

  const auto grid_shape = std::array { lattice.grid_extents().extent(0),
                                       lattice.grid_extents().extent(1),
                                       lattice.grid_extents().extent(2) };

  auto Ex = py::array_t<double, py::array::c_style>(grid_shape);
  auto Ey = py::array_t<double, py::array::c_style>(grid_shape);
  auto Ez = py::array_t<double, py::array::c_style>(grid_shape);
  auto Bx = py::array_t<double, py::array::c_style>(grid_shape);
  auto By = py::array_t<double, py::array::c_style>(grid_shape);
  auto Bz = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jx = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jy = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jz = py::array_t<double, py::array::c_style>(grid_shape);

  auto Exv = Ex.template mutable_unchecked<3>();
  auto Eyv = Ey.template mutable_unchecked<3>();
  auto Ezv = Ez.template mutable_unchecked<3>();
  auto Bxv = Bx.template mutable_unchecked<3>();
  auto Byv = By.template mutable_unchecked<3>();
  auto Bzv = Bz.template mutable_unchecked<3>();
  auto Jxv = Jx.template mutable_unchecked<3>();
  auto Jyv = Jy.template mutable_unchecked<3>();
  auto Jzv = Jz.template mutable_unchecked<3>();

  for(const auto mds = lattice.mds(); const auto idx: tyvi::sstd::index_space(mds)) {
    const auto [i, j, k] = idx;
    const auto F         = mds[idx][];

    Exv(i, j, k) = F.Ex;
    Eyv(i, j, k) = F.Ey;
    Ezv(i, j, k) = F.Ez;
    Bxv(i, j, k) = F.Bx;
    Byv(i, j, k) = F.By;
    Bzv(i, j, k) = F.Bz;
    Jxv(i, j, k) = F.Jx;
    Jyv(i, j, k) = F.Jy;
    Jzv(i, j, k) = F.Jz;
  }

  return std::tuple { std::tuple { std::move(Ex), std::move(Ey), std::move(Ez) },
                      std::tuple { std::move(Bx), std::move(By), std::move(Bz) },
                      std::tuple { std::move(Jx), std::move(Jy), std::move(Jz) } };
}
}  // namespace


//--------------------------------------------------

namespace emf {

void
  bind_emf(py::module& m_sub)
{

  //--------------------------------------------------
  // 1D bindings
  // TODO

  //--------------------------------------------------
  // 2D bindings
  // TODO

  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  using antenna_pyvec = pybind11::array_t<emf::antenna_mode::value_type>;
  using antenna_complex_pyvec =
    pybind11::array_t<std::complex<emf::antenna_mode::value_type>>;
  py::class_<emf::antenna_mode>(m_3d, "antenna_mode")
    .def(
      py::init([](
                 antenna_pyvec A,
                 std::optional<antenna_pyvec> k,
                 std::optional<antenna_pyvec> n,
                 std::optional<antenna_complex_pyvec> lap_coeffs) {
        if((k and n) or (not k and not n)) {
          throw std::runtime_error(
            "antenna_mode expects k or n to be defined but not both.");
        }


        auto assert_3d_vec = [](auto& x) {
          if(x.ndim() != 1) {
            throw std::runtime_error(
              "Antenna expects A and k/n to be rank-1 arrays (specifically 3D "
              "vectors).");
          }

          if(x.shape(0) != 3) {
            throw std::runtime_error("Antenna expects A and k/n to be 3D vectors.");
          }
        };

        const auto wave_data =
          std::invoke([&] -> decltype(emf::antenna_mode::wave_data) {
            if(k) {
              assert_3d_vec(k.value());
              const auto kv = k.value().template unchecked<1>();
              return emf::antenna_mode::wave_vector { { kv(0), kv(1), kv(2) } };
            } else {
              assert_3d_vec(n.value());
              const auto nv = n.value().template unchecked<1>();
              return emf::antenna_mode::wave_number { { nv(0), nv(1), nv(2) } };
            }
          });

        assert_3d_vec(A);
        const auto Av = A.template unchecked<1>();

        auto to_stdvec = [](antenna_complex_pyvec& p)
          -> std::optional<std::vector<std::complex<emf::antenna_mode::value_type>>> {
          if(p.ndim() != 1) {
            throw std::runtime_error { "lap_coeffs must be 1D array." };
          }

          const auto N  = static_cast<std::size_t>(p.shape(0));
          auto vec      = std::vector<std::complex<emf::antenna_mode::value_type>>(N);
          const auto pv = p.template unchecked<1>();

          for(auto i = 0uz; i < N; ++i) { vec[i] = pv(i); }
          return vec;
        };

        return emf::antenna_mode { .A { Av(0), Av(1), Av(2) },
                                   .wave_data { wave_data },
                                   .lap_coeffs { lap_coeffs.and_then(to_stdvec) } };
      }),
      py::kw_only(),
      py::arg("A"),
      py::arg("k")          = std::optional<antenna_pyvec> {},
      py::arg("n")          = std::optional<antenna_pyvec> {},
      py::arg("lap_coeffs") = std::optional<antenna_complex_pyvec> {});

  // 3d tile
  py::class_<emf::Tile<3>, corgi::Tile<3>, std::shared_ptr<emf::Tile<3>>>(m_3d, "Tile")
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return emf::Tile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def("set_EBJ", &emf::Tile<3>::set_EBJ)
    .def("batch_set_EBJ", &emf::Tile<3>::batch_set_EBJ)
    .def(
      "get_EBJ",
      [](emf::Tile<3>& tile) {
        auto EBJ = tile.get_EBJ();
        return to_ndarrays(EBJ);
      })
    .def(
      "get_EBJ_with_halo",
      [](emf::Tile<3>& tile) {
        auto EBJ = tile.get_EBJ_with_halo();
        return to_ndarrays(EBJ);
      })
    .def("push_half_b", &emf::Tile<3>::push_half_b)
    .def("push_e", &emf::Tile<3>::push_e)
    .def("filter_current", &emf::Tile<3>::filter_current)
    .def(
      "global_coordinate_map",
      [](const emf::Tile<3>& tile) {
        return std::function<std::array<double, 3>(std::array<double, 3>)>(
          [map = tile.global_coordinate_map()](const std::array<double, 3> idx) {
            return map(idx[0], idx[1], idx[2]);
          });
      })
    .def("register_antenna", &emf::Tile<3>::register_antenna)
    .def("deposit_antenna_current", &emf::Tile<3>::deposit_antenna_current)
    .def("add_current", &emf::Tile<3>::add_current);


  // HDF5 snapshot writer
  py::class_<hdf5::FieldsWriter<3>>(m_3d, "Hdf5FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &hdf5::FieldsWriter<3>::write);

  // MPI-IO snapshot writer
  py::class_<mpiio::FieldsWriter<3>>(m_3d, "MpiioFieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int>(),
         py::arg("prefix"), py::arg("Nx"), py::arg("NxMesh"),
         py::arg("Ny"), py::arg("NyMesh"),
         py::arg("Nz"), py::arg("NzMesh"),
         py::arg("stride"), py::arg("nspecies") = 2)
    .def("write", &mpiio::FieldsWriter<3>::write)
    .def("write_collective", &mpiio::FieldsWriter<3>::write_collective);

  m_3d.def("_write_average_B_energy_density", &emf::write_average_B_energy_density)
    .def("_write_average_E_energy_density", &emf::write_average_E_energy_density);
}
}  // namespace emf
