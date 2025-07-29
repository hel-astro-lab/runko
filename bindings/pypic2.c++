#include "core/particles_common.h"
#include "core/pic2/tile.h"
#include "io/tasker.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"

#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace {
namespace py = pybind11;
template<typename T>
auto
  to_ndarray(const std::vector<T>& vec)
{
  const auto grid_shape = std::array { vec.size() };
  auto mda              = py::array_t<T, py::array::c_style>(grid_shape);
  auto mda_mut          = mda.template mutable_unchecked<1>();

  for(const auto i: std::views::iota(0uz, vec.size())) { mda_mut(i) = vec[i]; }

  return std::move(mda);
}
}  // namespace

//--------------------------------------------------
// experimental PIC v2 module


namespace pic2 {
namespace py = pybind11;

// python bindings for plasma classes & functions
void
  bind_pic2(py::module& m_sub)
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

  py::class_<
    pic2::Tile<3>,
    emf2::Tile<3>,
    corgi::Tile<3>,
    std::shared_ptr<pic2::Tile<3>>>(m_3d, "Tile")
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return pic2::Tile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def(
      "get_positions",
      [](pic2::Tile<3>& tile, const std::size_t p) {
        const auto [x, y, z] = tile.get_positions(p);
        return std::tuple { to_ndarray(x), to_ndarray(y), to_ndarray(z) };
      })
    .def(
      "get_velocities",
      [](pic2::Tile<3>& tile, const std::size_t p) {
        const auto [x, y, z] = tile.get_velocities(p);
        return std::tuple { to_ndarray(x), to_ndarray(y), to_ndarray(z) };
      })
    .def("inject_to_each_cell", &pic2::Tile<3>::inject_to_each_cell)
    .def("inject", &pic2::Tile<3>::inject)
    .def("push_particles", &pic2::Tile<3>::push_particles)
    .def("deposit_current", &pic2::Tile<3>::deposit_current);

  //--------------------------------------------------
  // Full IO

  // 1D
  // TODO

  // 2D
  // TODO

  // 3D
  m_3d.def("write_particles", &pic2::write_particles<3>);
  m_3d.def("read_particles", &pic2::read_particles<3>);
}

}  // namespace pic2
