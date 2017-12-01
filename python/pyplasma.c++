#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// #include "../definitions.h"
#include "../solvers.h"
#include "../solvers/momentumLagrangianSolver.c++"
#include "../solvers/spatialLagrangianSolver.c++"
#include "../cell.h"
#include "../grid.h"
#include "../maxwell.h"


/// trampoline class for VlasovVelocitySolver
class PyVlasovVelocitySolver : public vlasov::VlasovVelocitySolver {
  public:
    using vlasov::VlasovVelocitySolver::setCell;
    using vlasov::VlasovVelocitySolver::setInterpolator;
    void solve() override {
      PYBIND11_OVERLOAD_PURE(void, vlasov::VlasovVelocitySolver, solve, );
    }
};


/// trampoline class for VlasovSpatialSolver
class PyVlasovSpatialSolver : public vlasov::VlasovSpatialSolver {
  public:
    using vlasov::VlasovSpatialSolver::setGrid;
    using vlasov::VlasovSpatialSolver::setTargetCell;
    void solve() override {
      PYBIND11_OVERLOAD_PURE(void, vlasov::VlasovSpatialSolver, solve, );
    }
};



// python bindings for plasma classes & functions
PYBIND11_MODULE(pyplasma, m) {


  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");
  py::class_<vlasov::VlasovCell, 
             corgi::Cell, 
             std::shared_ptr<vlasov::VlasovCell>
             >(m, "VlasovCell")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("addData",  &vlasov::VlasovCell::addData)
    .def("getData",  &vlasov::VlasovCell::getData)
    // .def("getDataPtr",  &vlasov::VlasovCell::getDataPtr) // TODO needs to be shared_ptr 
    .def("clip",     &vlasov::VlasovCell::clip)
    .def("bark",     &vlasov::VlasovCell::bark);


  // Loading node bindings from corgi library
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<vlasov::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    .def("cycle",    &vlasov::Grid::cycle)
    .def("howl",     &vlasov::Grid::howl);



  // trampoline base class followed by the actual solver implementations
  // Momentum dimension solvers
  py::class_<vlasov::VlasovVelocitySolver, PyVlasovVelocitySolver> vvsol(m, "VlasovVelocitySolver" );
  vvsol
    .def(py::init<>())
    .def("setCell",         &vlasov::VlasovVelocitySolver::setCell)
    .def("setInterpolator", &vlasov::VlasovVelocitySolver::setInterpolator)
    .def("solve",           &vlasov::VlasovVelocitySolver::solve);

  py::class_<vlasov::MomentumLagrangianSolver>(m, "MomentumLagrangianSolver", vvsol)
    .def(py::init<>());


  // trampoline base class followed by the actual solver implementations
  // Spatial dimension solvers
  py::class_<vlasov::VlasovSpatialSolver, PyVlasovSpatialSolver> vssol(m, "VlasovSpatialSolver" );
  vssol
    .def(py::init<>())
    .def("setTargetCell" ,  &vlasov::VlasovSpatialSolver::setTargetCell)
    .def("setGrid",         &vlasov::VlasovSpatialSolver::setGrid)
    .def("solve",           &vlasov::VlasovSpatialSolver::solve);

  py::class_<vlasov::SpatialLagrangianSolver2nd>(m, "SpatialLagrangianSolver2nd", vssol)
    .def(py::init<>());


  /// General class for handling Maxwell's equations
  py::class_<maxwell::PlasmaCell>(m, "PlasmaCell")
    .def(py::init<size_t, size_t, size_t>())
    .def("pushE",    &maxwell::PlasmaCell::pushE)
    .def("pushHalfB",    &maxwell::PlasmaCell::pushHalfB);




}



