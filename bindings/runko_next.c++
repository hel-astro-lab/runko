#include "core/communication_common.h"
#include "core/particles_common.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"


//--------------------------------------------------

namespace runko_next {

namespace py = pybind11;

void
  bind_runko_next([[maybe_unused]] py::module& m_sub)
{
  using PS = runko::ParticleState;
  py::class_<PS>(m_sub, "ParticleState")
    .def(py::init<PS::vec3, PS::vec3>(), py::arg("pos"), py::arg("vel"))
    .def_readwrite("pos", &PS::pos)
    .def_readwrite("vel", &PS::vel);

  m_sub.def("_virtual_tile_sync_handshake_mode", &runko::virtual_tile_sync_handshake_mode);
}
}  // namespace runko_next
