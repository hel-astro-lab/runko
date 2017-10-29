#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../sheets.h"
#include "../bundles.h"
#include "../velomesh.h"


using namespace sheets;
using namespace vmesh;


// python bindings for various helper classes 
PYBIND11_MODULE(plasmatools, m) {

  // --------------------------------------------------
  // Python bindings for sheet class
  py::class_<Sheet>(m, "Sheet" )
    .def(py::init<>())
    .def_readwrite("iGrid",      &Sheet::iGrid)
    .def_readwrite("jGrid",      &Sheet::jGrid)
    .def_readwrite("Ni",         &Sheet::Ni)
    .def_readwrite("Nj",         &Sheet::Nj)
    .def_readwrite("values",     &Sheet::values)
    .def("__iadd__",             &Sheet::operator+= )
    .def("__add__",              [](Sheet const & self, Sheet const other) 
        { return self + other; }, py::is_operator())
  .def("__isub__",             &Sheet::operator-= )
    .def("__sub__",              [](Sheet const & self, Sheet const other) 
        { return self - other; }, py::is_operator())
  .def("__imul__",             &Sheet::operator*= )
    .def("__mul__",              [](Sheet const & self, Realf const other) 
        { return self * other; }, py::is_operator())
  .def("resize",               &Sheet::resize)
    .def("loadValue",            &Sheet::loadValue)
    .def("getBlock",             &Sheet::getBlock)
    .def("isNonZero",            &Sheet::isNonZero);


  // --------------------------------------------------
  // Bundle bindings
  // TODO expand
  py::class_<bundles::Bundle>(m, "Bundle" )
    .def(py::init<>())
    .def("getGrid",   &bundles::Bundle::getGrid)
    .def("getPencil", &bundles::Bundle::getPencil);



  // --------------------------------------------------
  // VeloMesh bindings
  py::class_<vmesh::VeloMesh>(m, "VeloMesh" )
    .def(py::init<>())
    .def_readwrite("number_of_blocks", &vmesh::VeloMesh::number_of_blocks)
    .def_readwrite("mins",             &vmesh::VeloMesh::mins)
    .def_readwrite("maxs",             &vmesh::VeloMesh::maxs)
    .def_readwrite("lens",             &vmesh::VeloMesh::lens)
    .def_readwrite("Nblocks",          &vmesh::VeloMesh::Nblocks)
    .def("zFill",                      &vmesh::VeloMesh::zFill)
    .def("get_block",                  &vmesh::VeloMesh::get_block)
    .def("get_block_ID",               &vmesh::VeloMesh::get_block_ID)
    .def("get_indices",                &vmesh::VeloMesh::get_indices)
    .def("all_blocks",                 &vmesh::VeloMesh::all_blocks)
    .def("get_size",                   &vmesh::VeloMesh::get_size)
    .def("get_center",                 &vmesh::VeloMesh::get_center)
    .def("get_bundle",                 &vmesh::VeloMesh::get_bundle)
    .def("add_bundle",                 &vmesh::VeloMesh::add_bundle)
    .def("getSheet",                   &vmesh::VeloMesh::getSheet)

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
  .def("sizeInBytes", &vmesh::VeloMesh::sizeInBytes)
    .def("capacityInBytes", &vmesh::VeloMesh::capacityInBytes)
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


}

