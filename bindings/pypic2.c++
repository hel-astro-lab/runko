#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


//--------------------------------------------------
// experimental PIC v2 module


namespace pic2 {

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
}

}  // namespace pic2
