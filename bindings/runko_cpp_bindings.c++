#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "core/communication_common.h"
#include "runko_cpp_bindings.h"

//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(runko_cpp_bindings, m_base) {
  m_base.doc() = "Runko Python3 bindings";

  auto pycorgi = py::module::import("pycorgi");

  /// auxiliary tools
  py::module m_tools = m_base.def_submodule("tools", "auxiliary tools");
  tools::bind_tools(m_tools);

  /// emf
  py::module m_emf = m_base.def_submodule("emf", "emf module");
  emf::bind_emf(m_emf);

  /// pic
  py::module m_pic = m_base.def_submodule("pic", "pic module");
  pic::bind_pic(m_pic);

  //--------------------------------------------------
  // deprecated modules; TODO: add them back
  //
  /// vlv
  //py::module m_vlv = m_base.def_submodule("vlv", "Runko Vlasov module");
  //vlv::bind_vlv(m_vlv);
    
  /// qed
  //py::module m_qed = m_base.def_submodule("qed", "Runko QED module");
  //qed::bind_qed(m_qed);

  /// ffe
  //py::module m_ffe = m_base.def_submodule("ffe", "Runko force-free MHD module");
  //ffe::bind_ffe(m_ffe);

}

