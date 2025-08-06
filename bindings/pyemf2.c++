#include "core/emf2/tile.h"
#include "io/snapshots/fields2.h"
#include "io/tasker.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <memory>
#include <tuple>

namespace {

namespace py = pybind11;

/// MVP implementation.
auto
  to_ndarrays(emf2::YeeLattice::YeeLatticeHostCopy& lattice)
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

namespace emf2 {

void
  bind_emf2(py::module& m_sub)
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

  py::class_<emf2::Tile<3>, corgi::Tile<3>, std::shared_ptr<emf2::Tile<3>>>(
    m_3d,
    "Tile")
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return emf2::Tile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def("set_EBJ", &emf2::Tile<3>::set_EBJ)
    .def("batch_set_EBJ", &emf2::Tile<3>::batch_set_EBJ)
    .def(
      "get_EBJ",
      [](emf2::Tile<3>& tile) {
        auto EBJ = tile.get_EBJ();
        return to_ndarrays(EBJ);
      })
    .def(
      "get_EBJ_with_halo",
      [](emf2::Tile<3>& tile) {
        auto EBJ = tile.get_EBJ_with_halo();
        return to_ndarrays(EBJ);
      })
    .def("push_half_b", &emf2::Tile<3>::push_half_b)
    .def("push_e", &emf2::Tile<3>::push_e)
    .def("filter_current", &emf2::Tile<3>::filter_current)
    .def("subtract_J_from_E", &emf2::Tile<3>::subtract_J_from_E);

  py::class_<h5io::FieldsWriter2<3>>(m_3d, "FieldsWriter2")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::FieldsWriter2<3>::write);

  //--------------------------------------------------
  // Full IO

  // 1D
  // TODO

  // 2D
  // TODO

  // 3D
  m_3d.def("write_grids", &emf2::write_grids<3>);
  m_3d.def("read_grids", &emf2::read_grids<3>);
}
}  // namespace emf2
