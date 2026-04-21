#include "runko_cpp_bindings.h"

#include "pybind11/operators.h"
#include "core/communication_common.h"
#include "core/particles_common.h"
#include "tools/gpu_memory.h"

#include <exception>

namespace tools{

void bind_tools(pybind11::module& m)
{
  // mpi mode
  m.def("_virtual_tile_sync_handshake_mode", &runko::virtual_tile_sync_handshake_mode);

  py::module m_3d = m.def_submodule("threeD", "3D specializations");

  // Common communication enums.
  py::enum_<runko::comm_mode>(m, "comm_mode")
     .value("emf_E", runko::comm_mode::emf_E)
     .value("emf_B", runko::comm_mode::emf_B)
     .value("emf_J", runko::comm_mode::emf_J)
     .value("pic_particle", runko::comm_mode::pic_particle)
     .value("pic_particle_extra", runko::comm_mode::pic_particle_extra)
     .value("emf_J_exchange", runko::comm_mode::emf_J_exchange);

  // Common particle enums.
  py::enum_<runko::particle>(m, "particle")
     .value("electron", runko::particle::electron)
     .value("ion", runko::particle::ion)
     .value("photon", runko::particle::photon);

  m.def("_get_gpu_mem_kB", &runko::get_gpu_mem_kB,
        "Return GPU device memory in use (kB), or -1 on CPU backend.");

}

} // end of namespace tools
