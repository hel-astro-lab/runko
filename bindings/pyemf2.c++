#include "core/emf2/tile.h"
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
    }));
}
}  // namespace emf2
