#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../sheets.h"

using namespace sheets;


// python bindings for various helper classes 
PYBIND11_MODULE(plasmatools, m) {

    // Python bindings for sheet class
    py::class_<Sheet>(m, "Sheet" )
        .def(py::init<>())
        .def_readwrite("iGrid",      &Sheet::iGrid)
        .def_readwrite("jGrid",      &Sheet::jGrid)
        .def_readwrite("Ni",         &Sheet::Ni)
        .def_readwrite("Nj",         &Sheet::Nj)
        .def_readwrite("values",     &Sheet::values)
        .def("__iadd__",             &Sheet::operator+= )
        .def("__add__",           [](Sheet const & self, Sheet const other) 
            { return self + other; }, py::is_operator())
        .def("resize",               &Sheet::resize)
        .def("loadValue",            &Sheet::loadValue)
        .def("getBlock",             &Sheet::getBlock)
        .def("isNonZero",            &Sheet::isNonZero);



    // TODO add sheet arithmetics (+, -, /, *,...)
    
    /*
     * clsExampleOne.def("__eq__", &ExampleOne::operator==, py::is_operator());
     * clsExampleOne.def("__ne__", &ExampleOne::operator!=, py::is_operator());
     * clsExampleOne.def("__iadd__", &ExampleOne::operator+= );
     * clsExampleOne.def("__add__", [](ExampleOne const & self, ExampleOne const & other) { return self + other; }, py::is_operator());
     */
     
    
    
    
}

