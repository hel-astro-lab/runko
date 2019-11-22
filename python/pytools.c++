#include "py_submodules.h"
#include <pybind11/operators.h>

#include "../definitions.h"
#include "../tools/mesh.h"
#include "../vlasov/amr/mesh.h"

//#include "../em-fields/filters/filters.h"

#include "../tools/hilbert.h"

#include <exception>


namespace tools{

//--------------------------------------------------
// Different solver orders
using AM1d = toolbox::AdaptiveMesh<Realf, 1>;
using AM3d = toolbox::AdaptiveMesh<Realf, 3>;




// generator for Mesh bindings with type T and halo H
template<typename T, int H>
void declare_mesh(
    py::module &m, 
    const std::string& pyclass_name) 
{

    using Class = toolbox::Mesh<T, H>;
    //py::class_<Class>(m, pyclass_name.c_str())

    py::class_<
      toolbox::Mesh<T,H>,
      //std::shared_ptr<toolbox::Mesh<T,H>>,
      std::unique_ptr<toolbox::Mesh<T,H>,py::nodelete>
            >(m, pyclass_name.c_str())
    .def(py::init<size_t, size_t, size_t>())
    //.def("Nx", &Class::Nx)
    //.def("Ny", &Class::Ny)
    //.def("Nz", &Class::Nz)
    .def_property("Nx", [](Class &s){return (int)s.Nx;}, [](Class &s, int v){s.Nx = (size_t)v;})
    .def_property("Ny", [](Class &s){return (int)s.Ny;}, [](Class &s, int v){s.Ny = (size_t)v;})
    .def_property("Nz", [](Class &s){return (int)s.Nz;}, [](Class &s, int v){s.Nz = (size_t)v;})
    .def("indx",         &Class::indx)
    .def("size",         &Class::size)
    .def("__getitem__", [](Class &s, py::tuple indx) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();
        std::cout << "__getitem__ " << i << " " << j << " " << k << "\n";

        // NOTE: these are out-of-bounds; not inbound checks
        try {
        if (i < -H) throw py::index_error();
        if (j < -H) throw py::index_error();
        if (k < -H) throw py::index_error();

        if (i >= (int)s.Nx+H) throw py::index_error();
        if (j >= (int)s.Ny+H) throw py::index_error();
        if (k >= (int)s.Nz+H) throw py::index_error();
        } catch (std::exception& e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
        }
        std::cout << "__getitem passed ind__ nxnynz:" << s.Nx+H << " " << s.Ny+H << " " << s.Nz+H << "\n";

        T val = s(i,j,k);
        return val;

        return s(i,j,k);
      }, py::return_value_policy::reference)
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
    .def("clear",        &Class::clear)
    .def(py::self +  py::self)
    .def(py::self += py::self)
    .def(py::self -  py::self)
    .def(py::self -= py::self);
    //.def(py::self *  py::self)
    //.def(py::self *= py::self)
    //.def(py::self /  py::self)
    //.def(py::self /= py::self);
}




void bind_tools(pybind11::module& m)
{

  // declare Mesh with various halo sizes
  declare_mesh<Realf, 0>(m, "Mesh_H0" );
  declare_mesh<Realf, 1>(m, "Mesh_H1" );
  declare_mesh<Realf, 3>(m, "Mesh_H3" );



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
    .def("clip_neighbors",          &AM3d::clip_neighbors)
    .def("is_leaf",                 &AM3d::is_leaf)
    .def("set_min",                 &AM3d::set_min)
    .def("set_max",                 &AM3d::set_max)
    .def("get_size",                &AM3d::get_size)
    .def("get_length",              &AM3d::get_length)
    .def("get_center",              &AM3d::get_center)
    .def("get_level_0_cell_length", &AM3d::get_level_0_cell_length);





  //--------------------------------------------------
  // 2D bindings

  py::module m_2d = m.def_submodule("twoD", "2D specializations");

  //py::class_<fields::Filter>(m_2d, "Filter")
  //  .def(py::init<int, int>())
  //  .def("init_kernel",             &fields::Filter::init_kernel)
  //  .def("init_gaussian_kernel",    &fields::Filter::init_gaussian_kernel)
  //  //.def("init_sinc_kernel",      &fields::Filter::init_sinc_kernel)
  //  .def("init_lowpass_fft_kernel", &fields::Filter::init_lowpass_fft_kernel)
  //  .def("init_3point",             &fields::Filter::init_3point_kernel)
  //  .def("fft_kernel",              &fields::Filter::fft_kernel)
  //  .def("fft_image_forward",       &fields::Filter::fft_image_forward)
  //  .def("fft_image_backward",      &fields::Filter::fft_image_backward)
  //  .def("apply_kernel",            &fields::Filter::apply_kernel)
  //  .def("get_padded_current",      &fields::Filter::get_padded_current)
  //  .def("set_current",             &fields::Filter::set_current)
  //  .def("direct_convolve_3point",  &fields::Filter::direct_convolve_3point)
  //  .def("set_image",               &fields::Filter::set_image)
  //  .def("set_kernel",              &fields::Filter::set_kernel)
  //  .def("get_kernel",              &fields::Filter::get_kernel, py::return_value_policy::reference)
  //  .def("get_image",               &fields::Filter::get_image,  py::return_value_policy::reference);




  py::class_<hilbert::Hilbert2D>(m_2d, "HilbertGen")
    .def(py::init<int, int>())
    .def("hindex", &hilbert::Hilbert2D::hindex)
    .def("inv",    &hilbert::Hilbert2D::inv);


  // TODO: 3D Hilbert generator


}

} // end of namespace tools
