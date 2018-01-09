#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "mesh.h"



PYBIND11_MODULE(pyplasmaDev, m) {

  typedef float Realf;
  // typedef toolbox::AdaptiveMesh<Realf, 1> AM1d;
  typedef toolbox::AdaptiveMesh<Realf, 3> AM3d;



  py::class_<AM3d >(m, "AdaptiveMesh3D")
    .def(py::init<>())
    .def_readwrite("length", &AM3d::length)
    .def("resize",  &AM3d::resize)
    .def("__getitem__", [](const AM3d &s, py::tuple indx) 
        { 
        int rfl    = indx[0].cast<int>();
        uint64_t i = indx[1].cast<uint64_t>();
        uint64_t j = indx[2].cast<uint64_t>();
        uint64_t k = indx[3].cast<uint64_t>();
        uint64_t cid = s.get_cell({{i,j,k}}, rfl);

        return s.get(cid);
        })
    .def("__setitem__", [](AM3d &s, py::tuple indx, const Realf v) 
        { 
        int rfl    = indx[0].cast<int>();
        uint64_t i = indx[1].cast<uint64_t>();
        uint64_t j = indx[2].cast<uint64_t>();
        uint64_t k = indx[3].cast<uint64_t>();
        uint64_t cid = s.get_cell({{i,j,k}}, rfl);

        s.set(cid, v);
        });





};






