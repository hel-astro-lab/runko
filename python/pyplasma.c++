#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// #include "../definitions.h"
#include "../solvers.h"
#include "../solvers/SplittedLagrangian.c++"
#include "../cell.h"
#include "../grid.h"



/// trampoline class for VlasovVelocitySolver
class PyVlasovVelocitySolver : public vlasov::VlasovVelocitySolver {
  public:
    using vlasov::VlasovVelocitySolver::setMesh;
    using vlasov::VlasovVelocitySolver::setInterpolator;
    void solve() override {
      PYBIND11_OVERLOAD_PURE(void, vlasov::VlasovVelocitySolver, solve, );
    }
};





// python bindings for plasma classes & functions
PYBIND11_MODULE(pyplasma, m) {


  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");
  py::class_<vlasov::VCell>(m, "VCell", corgiCell)
    .def(py::init<size_t, size_t, int >())
    .def("addData",  &vlasov::VCell::addData)
    .def("getData",  &vlasov::VCell::getData)
    .def("bark",     &vlasov::VCell::bark);


  // Loading node bindings from corgi library
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<vlasov::Grid>(m, "Grid", corgiNode)
    .def(py::init<>())
    .def("cycle",    &vlasov::Grid::cycle)
    .def("howl",     &vlasov::Grid::howl);



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



