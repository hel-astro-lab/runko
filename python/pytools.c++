#include "py_submodules.h"


#include "../definitions.h"
#include "../tools/mesh.h"
#include "../vlasov/amr/mesh.h"


//--------------------------------------------------
// Different solver orders
using AM1d = toolbox::AdaptiveMesh<Realf, 1>;
using AM3d = toolbox::AdaptiveMesh<Realf, 3>;




// generator for Mesh bindings with type T and halo H
template<typename T, int H>
void declare_Mesh(
    py::module &m, 
    const std::string& pyclass_name) 
{
    using Class = toolbox::Mesh<T, H>;
    py::class_<Class>(m, pyclass_name.c_str())

    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &Class::Nx)
    .def_readwrite("Ny", &Class::Ny)
    .def_readwrite("Nz", &Class::Nz)
    .def("indx",         &Class::indx)
    .def("__getitem__", [](const Class &s, py::tuple indx) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();


        // NOTE: these are out-of-bounds; not inbound checks
        if (i < -H) throw py::index_error();
        if (j < -H) throw py::index_error();
        if (k < -H) throw py::index_error();

        if (i >= (int)s.Nx+H) throw py::index_error();
        if (j >= (int)s.Ny+H) throw py::index_error();
        if (k >= (int)s.Nz+H) throw py::index_error();

        return s(i,j,k);
      })
    .def("__setitem__", [](Class &s, py::tuple indx, Realf val) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();

        if (i < -H) throw py::index_error();
        if (j < -H) throw py::index_error();
        if (k < -H) throw py::index_error();

        if (i >= (int)s.Nx+H) throw py::index_error();
        if (j >= (int)s.Ny+H) throw py::index_error();
        if (k >= (int)s.Nz+H) throw py::index_error();

        s(i,j,k) = val;
        })
    .def("clear",        &Class::clear);
}




void bind_tools(pybind11::module& m)
{

  // declare Mesh with various halo sizes
  declare_Mesh<Realf, 0>(m, std::string("Mesh0") );
  declare_Mesh<Realf, 1>(m, std::string("Mesh1") );
  declare_Mesh<Realf, 3>(m, std::string("Mesh3") );



  //--------------------------------------------------

  py::class_<AM3d >(m, "AdaptiveMesh3D")
    .def(py::init<>())
    .def_readwrite("length",                     &AM3d::length)
    .def_readwrite("maximum_refinement_level",   &AM3d::maximum_refinement_level)
    .def_readwrite("top_refinement_level",       &AM3d::top_refinement_level)

    .def("resize",                &AM3d::resize)
    .def("get_cell_from_indices", &AM3d::get_cell_from_indices)
    .def("get_indices",           &AM3d::get_indices)
    .def("get_refinement_level",  &AM3d::get_refinement_level)
    .def("get_parent_indices",    &AM3d::get_parent_indices)
    .def("get_parent",            &AM3d::get_parent)

    .def("get_maximum_possible_refinement_level",&AM3d::get_maximum_possible_refinement_level)
    .def("set_maximum_refinement_level",         &AM3d::set_maximum_refinement_level)
    .def("get_level_0_parent_indices",           &AM3d::get_level_0_parent_indices)
    .def("get_level_0_parent",                   &AM3d::get_level_0_parent)
    .def("get_children",                         &AM3d::get_children)
    .def("get_siblings",                         &AM3d::get_siblings)
    .def("get_cells",                            &AM3d::get_cells)
    .def("__getitem__", [](const AM3d &s, py::tuple indx) 
        { 
        auto i = indx[0].cast<uint64_t>();
        auto j = indx[1].cast<uint64_t>();
        auto k = indx[2].cast<uint64_t>();
        auto    rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        return s.get_from_roots(cid);
        })
    .def("__setitem__", [](AM3d &s, py::tuple indx, Realf v) 
        { 
        auto i = indx[0].cast<uint64_t>();
        auto j = indx[1].cast<uint64_t>();
        auto k = indx[2].cast<uint64_t>();
        auto   rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        s.set(cid, v);
        })

    .def("clip_cells",              &AM3d::clip_cells)
    .def("is_leaf",                 &AM3d::is_leaf)
    .def("set_min",                 &AM3d::set_min)
    .def("set_max",                 &AM3d::set_max)
    .def("get_size",                &AM3d::get_size)
    .def("get_length",              &AM3d::get_length)
    .def("get_center",              &AM3d::get_center)
    .def("get_level_0_cell_length", &AM3d::get_level_0_cell_length);





}








