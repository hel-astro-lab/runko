#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../tools/bundles.h"
#include "../vlv/velomesh.h"
#include "../tools/mesh.h"



// trampoline class for BundleInterpolator (because of inheritance from Base class)
class PyBundleInterpolator : public toolbox::BundleInterpolator {
  public:
    using toolbox::BundleInterpolator::BundleInterpolator;
    using toolbox::BundleInterpolator::set_bundle;
    using toolbox::BundleInterpolator::get_bundle;
    toolbox::Bundle interpolate() override {
      PYBIND11_OVERLOAD_PURE(toolbox::Bundle, toolbox::BundleInterpolator, interpolate, );
    }
};



// python bindings for various helper classes 
PYBIND11_MODULE(plasmatools, m) {

  // --------------------------------------------------
  // Python bindings for sheet class
  /* TODO There is a template magick hiccup with these now
  py::class_<sheets::Sheet>(m, "Sheet" )
    .def(py::init<>())
    .def_readwrite("iGrid",      &sheets::Sheet::iGrid)
    .def_readwrite("jGrid",      &sheets::Sheet::jGrid)
    .def_readwrite("Ni",         &sheets::Sheet::Ni)
    .def_readwrite("Nj",         &sheets::Sheet::Nj)
    .def_readwrite("values",     &sheets::Sheet::values)
    .def("__iadd__",             &sheets::Sheet::operator+= )
    .def("__add__",              [](sheets::Sheet const & self, sheets::Sheet const other) 
        { return self + other; }, py::is_operator())
  .def("__isub__",               &sheets::Sheet::operator-= )
  .def("__sub__",              [](sheets::Sheet const & self, sheets::Sheet const other) 
        { return self - other; }, py::is_operator())
  // .def("__imul__",               &sheets::Sheet::operator*= )
  .def("__imul__",             [](sheets::Sheet const& self, sheets::Sheet const other)
        { return self * other; }, py::is_operator())
  .def("__mul__",              [](sheets::Sheet const & self, float_m const other) 
        { return self * other; }, py::is_operator())
  .def("resize",                 &sheets::Sheet::resize)
    .def("load_value",            &sheets::Sheet::load_value)
    .def("get_block",             &sheets::Sheet::get_block)
    .def("is_non_zero",            &sheets::Sheet::is_non_zero);
  */



  // --------------------------------------------------
  // Bundle bindings
  // TODO expand
  py::class_<toolbox::Bundle>(m, "Bundle" )
    .def(py::init<>())
    .def("get_grid",   &toolbox::Bundle::get_grid)
    .def("get_pencil", &toolbox::Bundle::get_pencil);


  py::class_<toolbox::BundleInterpolator, PyBundleInterpolator> bintp(m, "BundleInterpolator" );
  bintp
    .def(py::init<>())
    .def("set_bundle",   &toolbox::BundleInterpolator::set_bundle)
    .def("get_bundle",   &toolbox::BundleInterpolator::get_bundle)
    .def("interpolate", &toolbox::BundleInterpolator::interpolate);

  py::class_<toolbox::BundleInterpolator2nd>(m, "BundleInterpolator2nd", bintp)
    .def(py::init<>());

  py::class_<toolbox::BundleInterpolator4th>(m, "BundleInterpolator4th", bintp)
    .def(py::init<>());

  py::class_<toolbox::BundleInterpolator4PIC>(m, "BundleInterpolator4PIC", bintp)
    .def(py::init<>());


  // --------------------------------------------------
  // VeloMesh bindings
  py::class_<vmesh::VeloMesh>(m, "VeloMesh" )
    .def(py::init<>())
    .def_readwrite("number_of_blocks", &vmesh::VeloMesh::number_of_blocks)
    .def_readwrite("mins",             &vmesh::VeloMesh::mins)
    .def_readwrite("maxs",             &vmesh::VeloMesh::maxs)
    .def_readwrite("lens",             &vmesh::VeloMesh::lens)
    .def_readwrite("Nblocks",          &vmesh::VeloMesh::Nblocks)
    .def("z_fill",                      &vmesh::VeloMesh::z_fill)
    .def("get_block",                   &vmesh::VeloMesh::get_block)
    .def("get_block_id",                 &vmesh::VeloMesh::get_block_id)
    .def("get_indices",                 &vmesh::VeloMesh::get_indices)
    .def("all_blocks",                  &vmesh::VeloMesh::all_blocks)
    .def("get_size",                    &vmesh::VeloMesh::get_size)
    .def("get_center",                  &vmesh::VeloMesh::get_center)
    .def("get_bundle",                  &vmesh::VeloMesh::get_bundle)
    .def("add_bundle",                  &vmesh::VeloMesh::add_bundle)
    .def("get_sheet",                   &vmesh::VeloMesh::get_sheet)

    .def("__getitem__", [](const vmesh::VeloMesh &s, uint64_t i) {
        return s.__getitem__(i);
        })
  .def("__setitem__", [](vmesh::VeloMesh &s, uint64_t i, vblock_t v) {
      s.__setitem__(i, v);
      })

  // i,j,k indexing based interface
  .def("__getitem__", [](const vmesh::VeloMesh &s, py::tuple indx) {
      size_t i = indx[0].cast<size_t>();
      size_t j = indx[1].cast<size_t>();
      size_t k = indx[2].cast<size_t>();
      return s.__getitem2__( i,j,k );
      })
  .def("__setitem__", [](vmesh::VeloMesh &s, py::tuple indx, vblock_t v) {
      size_t i = indx[0].cast<size_t>();
      size_t j = indx[1].cast<size_t>();
      size_t k = indx[2].cast<size_t>();
      return s.__setitem2__( i,j,k, v);
      })
  .def("size_in_bytes", &vmesh::VeloMesh::size_in_bytes)
    .def("capacity_in_bytes", &vmesh::VeloMesh::capacity_in_bytes)
    .def("clip", &vmesh::VeloMesh::clip);


  // Example of bare bones array interface
  /*
     .def("__getitem__", [](const Sequence &s, size_t i) {
     if (i >= s.size()) throw py::index_error();
     return s[i];
     })
     .def("__setitem__", [](Sequence &s, size_t i, float v) {
     if (i >= s.size()) throw py::index_error();
     s[i] = v;
     })
  // Slices [optional]
  .def("__getitem__", [](const Sequence &s, py::slice slice) -> Sequence* {
  size_t start, stop, step, slicelength;
  if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
  throw py::error_already_set();
  Sequence *seq = new Sequence(slicelength);
  for (size_t i = 0; i < slicelength; ++i) {
  (*seq)[i] = s[start]; start += step;
  }
  return seq;
  })
  .def("__setitem__", [](Sequence &s, py::slice slice, const Sequence &value) {
  size_t start, stop, step, slicelength;
  if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
  throw py::error_already_set();
  if (slicelength != value.size())
  throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
  for (size_t i = 0; i < slicelength; ++i) {
  s[start] = value[i]; start += step;
  }
  })        
  */


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


  py::class_<toolbox::Mesh<vmesh::VeloMesh,0>>(m, "VMesh")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &toolbox::Mesh<vmesh::VeloMesh,0>::Nx)
    .def_readwrite("Ny", &toolbox::Mesh<vmesh::VeloMesh,0>::Ny)
    .def_readwrite("Nz", &toolbox::Mesh<vmesh::VeloMesh,0>::Nz)
    .def("indx",         &toolbox::Mesh<vmesh::VeloMesh,0>::indx)
    .def("__getitem__", [](const toolbox::Mesh<vmesh::VeloMesh,0> &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        return s(i,j,k);
      }, py::return_value_policy::reference)
    .def("__setitem__", [](toolbox::Mesh<vmesh::VeloMesh,0> &s, py::tuple indx, vmesh::VeloMesh val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        s(i,j,k) = val;
        })
    .def("clear",        &toolbox::Mesh<vmesh::VeloMesh,0>::clear);

}

