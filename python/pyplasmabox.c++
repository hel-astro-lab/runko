#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "py_submodules.h"




//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(pyplasmabox, m_base) {
  m_base.doc() = "pyPlasmaBox Python3 bindings";

  auto pycorgi = py::module::import("pycorgi");
  //py::object corgiNode = (py::object) py::module::import("pycorgi").attr("Node");
  //py::object corgiTile = (py::object) py::module::import("pycorgi").attr("Tile");


  /// auxiliary tools
  py::module m_tools = m_base.def_submodule("tools", "PlasmaBox auxiliary tools");
  tools::bind_tools(m_tools);

  /// fields
  py::module m_fields = m_base.def_submodule("fields", "PlasmaBox fields module");
  fields::bind_fields(m_fields);

  /// vlv
  py::module m_vlv = m_base.def_submodule("vlv", "PlasmaBox Vlasov module");
  vlv::bind_vlv(m_vlv);

  /// pic
  py::module m_pic = m_base.def_submodule("pic", "PlasmaBox PIC module");
  pic::bind_pic(m_pic);

}

