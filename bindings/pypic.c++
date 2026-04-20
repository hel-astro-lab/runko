#include "core/particles_common.h"
#include "core/pic/reflector_wall.h"
#include "core/pic/tile.h"
#include "core/pic/virtual_tile.h"
#include "io/pic_average_kinetic_energy.h"
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
// experimental PIC module


namespace pic {
namespace py = pybind11;

// python bindings for plasma classes & functions
void
  bind_pic(py::module& m_sub)
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

  // object for storing single particle data
  using PS = runko::ParticleState<double>;
  py::class_<PS>(m_3d, "ParticleState")
    .def(py::init<PS::vec3, PS::vec3>(), py::arg("pos"), py::arg("vel"))
    .def_readwrite("pos", &PS::pos)
    .def_readwrite("vel", &PS::vel);

  // object for storing a batch of multiple particle's data
  using PSB = pic::ParticleStateBatch;
  py::class_<PSB>(m_3d, "ParticleStateBatch")
    .def(
      py::init<PSB::container_type, PSB::container_type>(),
      py::arg("pos"),
      py::arg("vel"))
    .def_readwrite("pos", &PSB::pos)
    .def_readwrite("vel", &PSB::vel);

  // reflector wall data structure
  using RW = pic::reflector_wall;
  py::class_<RW>(m_3d, "reflector_wall")
    .def(
      py::init(
        [](RW::value_type walloc, RW::value_type betawall, RW::value_type gammawall) {
          return RW { walloc, betawall, gammawall };
        }),
      py::arg("walloc"),
      py::arg("betawall")  = RW::value_type { 0 },
      py::arg("gammawall") = RW::value_type { 1 })
    .def_readwrite("walloc", &RW::walloc)
    .def_readwrite("betawall", &RW::betawall)
    .def_readwrite("gammawall", &RW::gammawall);

  // 3d virtual tile specialization
  py::class_<
    pic::VirtualTile<3>,
    emf::VirtualTile<3>,
    corgi::Tile<3>,
    std::shared_ptr<pic::VirtualTile<3>>>(m_3d, "VirtualTile")
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return pic::VirtualTile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def_static("canonical_type", []() { return py::type::of<pic::Tile<3>>(); });

  // 3d pic tile
  py::class_<pic::Tile<3>, emf::Tile<3>, corgi::Tile<3>, std::shared_ptr<pic::Tile<3>>>(
    m_3d,
    "Tile")
    .def_static("canonical_type", []() { return py::type::of<pic::Tile<3>>(); })
    .def_static(
      "virtual_tile_specialization",
      []() { return py::type::of<pic::VirtualTile<3>>(); })
    .def(
      py::init([](const std::array<std::size_t, 3> tile_grid_idx, const py::handle& h) {
        return pic::Tile<3>(tile_grid_idx, toolbox::ConfigParser(h));
      }))
    .def(
      "get_positions",
      [](pic::Tile<3>& tile, const std::size_t p) {
        const auto [x, y, z] = tile.get_positions(p);
        return std::tuple { to_ndarray(x), to_ndarray(y), to_ndarray(z) };
      })
    .def(
      "get_velocities",
      [](pic::Tile<3>& tile, const std::size_t p) {
        const auto [x, y, z] = tile.get_velocities(p);
        return std::tuple { to_ndarray(x), to_ndarray(y), to_ndarray(z) };
      })
    .def(
      "get_ids",
      [](pic::Tile<3>& tile, const std::size_t p) {
        return to_ndarray(tile.get_ids(p));
      })
    .def("pack_outgoing_particles", &pic::Tile<3>::pack_outgoing_particles)
    .def("inject_to_each_cell", &pic::Tile<3>::inject_to_each_cell)
    .def("inject", &pic::Tile<3>::inject)
    .def("batch_inject_to_cells", &pic::Tile<3>::batch_inject_to_cells)
    .def("push_particles", &pic::Tile<3>::push_particles)
    .def("deposit_current", &pic::Tile<3>::deposit_current)
    .def("sort_particles", &pic::Tile<3>::sort_particles)
    .def("inject_dead_particles", &pic::Tile<3>::inject_dead_particles)
    .def("register_reflector_wall", &pic::Tile<3>::register_reflector_wall)
    .def("reflect_particles", &pic::Tile<3>::reflect_particles)
    .def("advance_reflector_walls", &pic::Tile<3>::advance_reflector_walls)
    .def("batch_inject_in_x_stripe", &pic::Tile<3>::batch_inject_in_x_stripe);


  m_3d.def("_write_average_kinetic_energy", &pic::write_average_kinetic_energy);
}

}  // namespace pic
