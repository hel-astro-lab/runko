#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../sheets.h"

using namespace sheets;


// python bindings for various helper classes 
PYBIND11_MODULE(plasmaTools, m) {

    // Python bindings for sheet class
    py::class_<Sheet>(m, "Sheet" )
        .def(py::init<>())
        .def_readwrite("iGrid",      &Sheet::iGrid)
        .def_readwrite("jGrid",      &Sheet::jGrid)
        .def_readwrite("Ni",         &Sheet::Ni)
        .def_readwrite("Nj",         &Sheet::Nj)
        .def_readwrite("values",     &Sheet::values);

    // TODO add sheet arithmetics (+, -, /, *,...)
    
    
    
    
}

