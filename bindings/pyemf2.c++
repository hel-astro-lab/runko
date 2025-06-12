#include "core/emf2/tile.h"
#include "io/tasker.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"

#include <memory>
#include <string>
#include <unordered_map>


//--------------------------------------------------

namespace emf2 {

namespace py = pybind11;

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
    .def(py::init([](const py::handle& h) {
      return emf2::Tile<3>(toolbox::ConfigParser(h));
    }))
    .def("set_fields", &emf2::Tile<3>::set_fields);

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
