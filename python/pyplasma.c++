#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../definitions.h"
#include "../vlasov/cell.h"
#include "../vlasov/grid.h"
#include "../em-fields/fields.h"




// python bindings for plasma classes & functions
PYBIND11_MODULE(pyplasma, m) {

  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");


  /// General class for handling Maxwell's equations
  py::class_<fields::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<fields::PlasmaCell>
            >(m, "PlasmaCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("yeeDt",  &fields::PlasmaCell::yeeDt)
    .def_readwrite("yeeDx",  &fields::PlasmaCell::yeeDx)
    .def_readwrite("yeeDy",  &fields::PlasmaCell::yeeDy)
    .def_readwrite("yeeDz",  &fields::PlasmaCell::yeeDz)
    .def("cycleYee",         &fields::PlasmaCell::cycleYee)
    .def("pushE",            &fields::PlasmaCell::pushE)
    .def("pushHalfB",        &fields::PlasmaCell::pushHalfB)
    .def("depositCurrent",   &fields::PlasmaCell::depositCurrent)
    .def("getYee",           &fields::PlasmaCell::getYee, py::return_value_policy::reference)
    .def("updateBoundaries", &fields::PlasmaCell::updateBoundaries);



  py::class_<vlasov::VlasovCell, 
             fields::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<vlasov::VlasovCell>
             >(m, "VlasovCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",     &vlasov::VlasovCell::dt)
    .def_readwrite("dx",     &vlasov::VlasovCell::dx)
    .def_readwrite("dy",     &vlasov::VlasovCell::dy)
    .def_readwrite("dz",     &vlasov::VlasovCell::dz)
    .def("getPlasmaSpecies", [](vlasov::VlasovCell& cell, size_t i, size_t s) 
        { return cell.steps.get(i).at(s); }, py::return_value_policy::reference)
    .def("insertInitialSpecies", [](vlasov::VlasovCell& c, 
                                  std::vector<vlasov::PlasmaBlock> species){
        c.steps.push_back(species);
        c.steps.push_back(species);
        })

    .def("clip",         &vlasov::VlasovCell::clip)
    .def("cycle",        &vlasov::VlasovCell::cycle);



  // Loading node bindings from corgi library
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<vlasov::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>());


  py::class_<fields::YeeLattice>(m, "YeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex", &fields::YeeLattice::ex)
    .def_readwrite("ey", &fields::YeeLattice::ey)
    .def_readwrite("ez", &fields::YeeLattice::ez)
    .def_readwrite("bx", &fields::YeeLattice::bx)
    .def_readwrite("by", &fields::YeeLattice::by)
    .def_readwrite("bz", &fields::YeeLattice::bz)
    .def_readwrite("jx", &fields::YeeLattice::jx)
    .def_readwrite("jy", &fields::YeeLattice::jy)
    .def_readwrite("jz", &fields::YeeLattice::jz)
    .def_readwrite("rh", &fields::YeeLattice::rh);


}



