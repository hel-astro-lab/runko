#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// #include "../definitions.h"
#include "../solvers.h"
#include "../solvers/SplittedLagrangian.c++"



/// trampoline class for VlasovVelocitySolver
class PyVlasovVelocitySolver : public vlasov::VlasovVelocitySolver {
  public:
    using vlasov::VlasovVelocitySolver::setMesh;
    using vlasov::VlasovVelocitySolver::setInterpolator;
    void solve() override {
      PYBIND11_OVERLOAD_PURE(void, vlasov::VlasovVelocitySolver, solve, );
    }
};



// python bindings for plasma classes & function
PYBIND11_MODULE(pyplasma, m) {

  // trampoline base class followed by the actual solver implementations
  py::class_<vlasov::VlasovVelocitySolver, PyVlasovVelocitySolver> vsol(m, "VlasovVelocitySolver" );
  vsol
    .def(py::init<>())
    .def("setMesh",         &vlasov::VlasovVelocitySolver::setMesh)
    .def("setInterpolator", &vlasov::VlasovVelocitySolver::setInterpolator)
    .def("solve",           &vlasov::VlasovVelocitySolver::solve);

  py::class_<vlasov::SplittedLagrangian>(m, "SplittedLagrangian", vsol)
    .def(py::init<>());



  // py::class_<vmesh::sSolver>(m, "sSolver" )
  //   .def(py::init<vmesh::Node&>())
  //   // .def_readwrite("node",   &vmesh::sSolver::node)
  //   // .def("setNode" ,         &vmesh::sSolver::setNode)
  //   .def("setTargetCell" ,   &vmesh::sSolver::setTargetCell)
  //   .def("solve",            &vmesh::sSolver::solve);
  // // .def("update",           &vmesh::sSolver::update);



}



