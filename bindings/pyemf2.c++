#include "py_submodules.h"

#include <string>


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
}
}  // namespace emf2
