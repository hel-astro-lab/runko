#include "core/emf2/tile.h"
#include "io/tasker.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <memory>

namespace {

namespace py = pybind11;

/// MVP implementation. Only supports 3D vector grids.
template<
  typename V,
  typename ElemExtents,
  typename ElemLP,
  typename GridExtents,
  typename GridLP>
  requires(
    ElemExtents::rank() == 1uz and ElemExtents::static_extent(0) == 3uz and
    GridExtents::rank() == 3uz)
auto
  tyvi_mdgrid_buffer_to_ndarray(
    const tyvi::mdgrid_buffer<V, ElemExtents, ElemLP, GridExtents, GridLP>& mdgbuff)
{

  const auto grid_shape = std::array { mdgbuff.grid_extents().extent(0),
                                       mdgbuff.grid_extents().extent(1),
                                       mdgbuff.grid_extents().extent(2) };
  constexpr auto b      = sizeof(typename V::value_type);
  const auto grid_strides =
    std::array { b * grid_shape[1] * grid_shape[2], b * grid_shape[2], b };

  auto x = py::array_t<typename V::value_type>(grid_shape, grid_strides);
  auto y = py::array_t<typename V::value_type>(grid_shape, grid_strides);
  auto z = py::array_t<typename V::value_type>(grid_shape, grid_strides);

  auto [xv, yv, zv] = std::tuple { x.template mutable_unchecked<3>(),
                                   y.template mutable_unchecked<3>(),
                                   z.template mutable_unchecked<3>() };

  for(const auto mds = mdgbuff.mds(); const auto idx: tyvi::sstd::index_space(mds)) {
    const auto [i, j, k] = idx;

    xv(i, j, k) = mds[idx][0];
    yv(i, j, k) = mds[idx][1];
    zv(i, j, k) = mds[idx][2];
  }

  return std::tuple { std::move(x), std::move(y), std::move(z) };
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
    .def("get_EBJ", [](emf2::Tile<3>& tile) {
      const auto [E, B, J] = tile.get_EBJ();

      return std::tuple { tyvi_mdgrid_buffer_to_ndarray(E),
                          tyvi_mdgrid_buffer_to_ndarray(B),
                          tyvi_mdgrid_buffer_to_ndarray(J) };
    });

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
