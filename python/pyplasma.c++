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



  /// General class for handling Maxwell's equations
  py::class_<maxwell::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<maxwell::PlasmaCell>
            >(m, "PlasmaCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("yeeDt",  &maxwell::PlasmaCell::yeeDt)
    .def_readwrite("yeeDx",  &maxwell::PlasmaCell::yeeDx)
    .def_readwrite("yeeDy",  &maxwell::PlasmaCell::yeeDy)
    .def_readwrite("yeeDz",  &maxwell::PlasmaCell::yeeDz)
    .def("cycleYee",         &maxwell::PlasmaCell::cycleYee)
    .def("pushE",            &maxwell::PlasmaCell::pushE)
    .def("pushHalfB",        &maxwell::PlasmaCell::pushHalfB)
    .def("depositCurrent",   &maxwell::PlasmaCell::depositCurrent)
    .def("getYee",           &maxwell::PlasmaCell::getYee,    py::return_value_policy::reference)
    .def("getNewYee",        &maxwell::PlasmaCell::getNewYee, py::return_value_policy::reference)
    .def("updateBoundaries", &maxwell::PlasmaCell::updateBoundaries);



  py::class_<vlasov::VlasovCell, 
             maxwell::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<vlasov::VlasovCell>
             >(m, "VlasovCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",     &vlasov::VlasovCell::dt)
    .def_readwrite("dx",     &vlasov::VlasovCell::dx)
    .def_readwrite("dy",     &vlasov::VlasovCell::dy)
    .def_readwrite("dz",     &vlasov::VlasovCell::dz)
    .def("getPlasmaGrid",    &vlasov::VlasovCell::getPlasmaGrid, py::return_value_policy::reference)
    .def("getNewPlasmaGrid", &vlasov::VlasovCell::getNewPlasmaGrid, py::return_value_policy::reference)
    .def("clip",        &vlasov::VlasovCell::clip)
    .def("cyclePlasma", &vlasov::VlasovCell::cycle)
    .def("bark",        &vlasov::VlasovCell::bark);



  // Loading node bindings from corgi library
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<vlasov::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    //.def("cycle",    &vlasov::Grid::cycle)
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

  py::class_<vlasov::SpatialLagrangianSolver4th>(m, "SpatialLagrangianSolver4th", vssol)
    .def(py::init<>());


  py::class_<maxwell::YeeLattice>(m, "YeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex", &maxwell::YeeLattice::ex)
    .def_readwrite("ey", &maxwell::YeeLattice::ey)
    .def_readwrite("ez", &maxwell::YeeLattice::ez)
    .def_readwrite("bx", &maxwell::YeeLattice::bx)
    .def_readwrite("by", &maxwell::YeeLattice::by)
    .def_readwrite("bz", &maxwell::YeeLattice::bz)
    .def_readwrite("jx", &maxwell::YeeLattice::jx)
    .def_readwrite("jy", &maxwell::YeeLattice::jy)
    .def_readwrite("jz", &maxwell::YeeLattice::jz)
    .def_readwrite("rh", &maxwell::YeeLattice::rh);


  // --------------------------------------------------


  py::class_<vlasov::VlasovFluid>(m, "VlasovFluid")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("electrons", &vlasov::VlasovFluid::electrons)
    .def_readwrite("positrons", &vlasov::VlasovFluid::positrons)
    .def_readwrite("qms",       &vlasov::VlasovFluid::qms);


}



