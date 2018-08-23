#include "py_submodules.h"


#include "../definitions.h"
#include "../tools/mesh.h"



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





}








