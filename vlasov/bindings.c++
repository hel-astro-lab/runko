#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "mesh.h"
#include "amr_numerics.h"
#include "amr_refiner.h"



PYBIND11_MODULE(pyplasmaDev, m) {

  typedef float Realf;
  // typedef toolbox::AdaptiveMesh<Realf, 1> AM1d;
    
  typedef toolbox::AdaptiveMesh<Realf, 3> AM3d;
  typedef toolbox::Adapter<Realf, 3> Adapter3d;



  py::class_<AM3d >(m, "AdaptiveMesh3D")
    .def(py::init<>())
    .def_readwrite("length",                     &AM3d::length)
    .def_readwrite("maximum_refinement_level",   &AM3d::maximum_refinement_level)
    .def("resize",               &AM3d::resize)
    .def("get_cell",             &AM3d::get_cell_from_indices)
    .def("get_indices",          &AM3d::get_indices)
    .def("get_refinement_level", &AM3d::get_refinement_level)
    .def("get_parent_indices",   &AM3d::get_parent_indices)
    .def("get_parent",           &AM3d::get_parent)
    .def("get_maximum_possible_refinement_level",&AM3d::get_maximum_possible_refinement_level)
    .def("get_level_0_parent_indices", &AM3d::get_level_0_parent_indices)
    .def("get_level_0_parent",   &AM3d::get_level_0_parent)
    .def("get_children",         &AM3d::get_children)
    .def("get_cells",            &AM3d::get_cells)
    .def("__getitem__", [](const AM3d &s, py::tuple indx) 
        { 
        uint64_t i = indx[0].cast<uint64_t>();
        uint64_t j = indx[1].cast<uint64_t>();
        uint64_t k = indx[2].cast<uint64_t>();
        int    rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        return s.get_from_roots(cid);
        })
    .def("__setitem__", [](AM3d &s, py::tuple indx, Realf v) 
        { 
        uint64_t i = indx[0].cast<uint64_t>();
        uint64_t j = indx[1].cast<uint64_t>();
        uint64_t k = indx[2].cast<uint64_t>();
        int   rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        s.set(cid, v);
        })
    .def("clip_cells",       &AM3d::clip_cells)
    //.def("cut_roots",        &AM3d::cut_roots)
    .def("split",                   &AM3d::split)
    .def("update_from_children",    &AM3d::update_from_children)
    .def("update_from_leafs",       &AM3d::update_from_leafs)
    .def("get_from_leafs",          &AM3d::get_from_leafs)
    .def("is_leaf",                 &AM3d::is_leaf)
    .def("set_min",                 &AM3d::set_min)
    .def("set_max",                 &AM3d::set_max)
    .def("get_length",              &AM3d::get_length)
    .def("get_center",              &AM3d::get_center)
    .def("get_level_0_cell_length", &AM3d::get_level_0_cell_length);



    // -------------------------------------------------- 
    // AMR numerics
      
  m.def("deriv", &toolbox::deriv<Realf, 3>);

  m.def("grad",  &toolbox::grad<Realf, 3>);



  py::class_<Adapter3d>(m, "Adapter")
    .def(py::init<>())
    .def_readwrite("tolerance", &Adapter3d::tolerance)
    .def_readwrite("cells_to_refine",   &Adapter3d::cells_to_refine)
    .def_readwrite("cells_to_unrefine", &Adapter3d::cells_to_unrefine)
    .def_readwrite("cells_created",     &Adapter3d::cells_created)
    .def("check",     &Adapter3d::check)
    .def("refine",    &Adapter3d::refine);



};






