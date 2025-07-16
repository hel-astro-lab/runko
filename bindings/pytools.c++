#include "py_submodules.h"

#include "pybind11/operators.h"
#include "definitions.h"
#include "tools/mesh.h"
#include "core/communication_common.h"
#include "core/particles_common.h"
#include "core/vlv/amr/mesh.h"
#include "tools/hilbert.h"

#include <exception>


namespace tools{

//--------------------------------------------------
// Different solver orders
using AM1d = toolbox::AdaptiveMesh<float, 1>;
using AM3d = toolbox::AdaptiveMesh<float, 3>;




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
      std::shared_ptr<toolbox::Mesh<T,H>>
      //std::unique_ptr<toolbox::Mesh<T,H>,py::nodelete>
            >(m, pyclass_name.c_str())
    .def(py::init<int, int, int>())
    //.def("Nx", &Class::Nx)
    //.def("Ny", &Class::Ny)
    //.def("Nz", &Class::Nz)
    .def_property("Nx", [](Class &s){return s.Nx;}, [](Class &s, int v){s.Nx = v;})
    .def_property("Ny", [](Class &s){return s.Ny;}, [](Class &s, int v){s.Ny = v;})
    .def_property("Nz", [](Class &s){return s.Nz;}, [](Class &s, int v){s.Nz = v;})
    .def("get_Ny", [](Class &s){ return s.Ny;})
    .def("indx",         &Class::indx)
    .def("size",         &Class::size)
    .def("__getitem__", [](Class &s, const py::tuple& indx) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();

        // NOTE: these are out-of-bounds; not inbound checks
        try {
          if (i < -H) throw py::index_error();
          if (j < -H) throw py::index_error();
          if (k < -H) throw py::index_error();

          if (i >= (int)s.Nx+H) throw py::index_error();
          if (j >= (int)s.Ny+H) throw py::index_error();
          if (k >= (int)s.Nz+H) throw py::index_error();
        } catch (std::exception& e) {
          std::cerr << "Standard exception: " << e.what() << std::endl;
        }

        T val = s(i,j,k);
        return val;

        return s(i,j,k);
      }) //, py::return_value_policy::reference)
    .def("__setitem__", [](Class &s, const py::tuple& indx, float val) 
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
  declare_mesh<float, 0>(m, "Mesh_H0" );
  declare_mesh<float, 1>(m, "Mesh_H1" );
  declare_mesh<float, 3>(m, "Mesh_H3" );


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
    .def("__setitem__", [](AM3d &s, py::tuple indx, float v) 
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

  //py::class_<emf::Filter>(m_2d, "Filter")
  //  .def(py::init<int, int>())
  //  .def("init_kernel",             &emf::Filter::init_kernel)
  //  .def("init_gaussian_kernel",    &emf::Filter::init_gaussian_kernel)
  //  //.def("init_sinc_kernel",      &emf::Filter::init_sinc_kernel)
  //  .def("init_lowpass_fft_kernel", &emf::Filter::init_lowpass_fft_kernel)
  //  .def("init_3point",             &emf::Filter::init_3point_kernel)
  //  .def("fft_kernel",              &emf::Filter::fft_kernel)
  //  .def("fft_image_forward",       &emf::Filter::fft_image_forward)
  //  .def("fft_image_backward",      &emf::Filter::fft_image_backward)
  //  .def("apply_kernel",            &emf::Filter::apply_kernel)
  //  .def("get_padded_current",      &emf::Filter::get_padded_current)
  //  .def("set_current",             &emf::Filter::set_current)
  //  .def("direct_convolve_3point",  &emf::Filter::direct_convolve_3point)
  //  .def("set_image",               &emf::Filter::set_image)
  //  .def("set_kernel",              &emf::Filter::set_kernel)
  //  .def("get_kernel",              &emf::Filter::get_kernel, py::return_value_policy::reference)
  //  .def("get_image",               &emf::Filter::get_image,  py::return_value_policy::reference);



  // 2D Hilbert generator
  py::class_<hilbert::Hilbert2D>(m_2d, "HilbertGen")
    .def(py::init<int, int>())
    .def("hindex", &hilbert::Hilbert2D::hindex)
    .def("inv",    &hilbert::Hilbert2D::inv);


  py::module m_3d = m.def_submodule("threeD", "3D specializations");

  // 3D Hilbert generator
  py::class_<hilbert::Hilbert3D>(m_3d, "HilbertGen")
    .def(py::init<int, int, int>())
    .def("hindex", &hilbert::Hilbert3D::hindex)
    .def("inv",    &hilbert::Hilbert3D::inv);

  // Common communication enums.
  py::enum_<runko::comm_mode>(m, "comm_mode")
     .value("emf_E", runko::comm_mode::emf_E)
     .value("emf_B", runko::comm_mode::emf_B)
     .value("emf_J", runko::comm_mode::emf_J)
     .value("pic_particle", runko::comm_mode::pic_particle)
     .value("pic_particle_extra", runko::comm_mode::pic_particle_extra);

  // Common particle enums.
  py::enum_<runko::particle>(m, "particle")
     .value("electron", runko::particle::electron)
     .value("ion", runko::particle::ion)
     .value("photon", runko::particle::photon);

}

} // end of namespace tools
