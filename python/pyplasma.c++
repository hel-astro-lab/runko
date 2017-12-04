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
    .def("pushE",            &maxwell::PlasmaCell::pushE)
    .def("pushHalfB",        &maxwell::PlasmaCell::pushHalfB)
    .def("getYee",           &maxwell::PlasmaCell::getYee,    py::return_value_policy::reference)
    .def("getNewYee",        &maxwell::PlasmaCell::getNewYee, py::return_value_policy::reference)
    .def("updateBoundaries", &maxwell::PlasmaCell::updateBoundaries);
  



  py::class_<vlasov::VlasovCell, 
             maxwell::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<vlasov::VlasovCell>
             >(m, "VlasovCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
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

  py::class_<toolbox::Mesh<double,1>>(m, "Mesh")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &toolbox::Mesh<double,1>::Nx)
    .def_readwrite("Ny", &toolbox::Mesh<double,1>::Ny)
    .def_readwrite("Nz", &toolbox::Mesh<double,1>::Nz)
    .def("indx",         &toolbox::Mesh<double,1>::indx)
    .def("__getitem__", [](const toolbox::Mesh<double,1> &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < -1) throw py::index_error();
        if (j < -1) throw py::index_error();
        if (k < -1) throw py::index_error();

        if (i > (int)s.Nx+1) throw py::index_error();
        if (j > (int)s.Ny+1) throw py::index_error();
        if (k > (int)s.Nz+1) throw py::index_error();

        return s(i,j,k);
      })
    .def("__setitem__", [](toolbox::Mesh<double,1> &s, py::tuple indx, double val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < -1) throw py::index_error();
        if (j < -1) throw py::index_error();
        if (k < -1) throw py::index_error();

        if (i > (int)s.Nx+1) throw py::index_error();
        if (j > (int)s.Ny+1) throw py::index_error();
        if (k > (int)s.Nz+1) throw py::index_error();

        s(i,j,k) = val;
        })
    .def("clear",        &toolbox::Mesh<double,1>::clear);

  py::class_<maxwell::YeeLattice>(m, "YeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    // .def("ex", &maxwell::YeeLattice::getEx, py::return_value_policy::reference)
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




}



