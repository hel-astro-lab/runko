#include "core/pic2/tile.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "tools/config_parser.h"

#include <memory>
#include <string>
#include <unordered_map>


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
  [[maybe_unused]] py::module m_3d =
    m_sub.def_submodule("threeD", "3D specializations");

  py::class_<
    pic2::Tile<3>,
    emf2::Tile<3>,
    corgi::Tile<3>,
    std::shared_ptr<pic2::Tile<3>>>(m_3d, "Tile", py::multiple_inheritance())
    .def(py::init([](const py::handle& h) {
      return pic2::Tile<3>(toolbox::ConfigParser(h));
    }));
}

}  // namespace pic2
