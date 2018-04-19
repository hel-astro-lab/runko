#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


//--------------------------------------------------
// experimental PIC module
  
#include "../pic/cell.h"
#include "../pic/pusher.h"
#include "../pic/field_interpolator.h"





//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(pypic, m) {

  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  py::object plasmaCell = (py::object) py::module::import("pyplasma").attr("PlasmaCell");




  py::class_<pic::PicCell, 
             fields::PlasmaCell,
             corgi::Cell, 
             std::shared_ptr<pic::PicCell>
             >(m, "PicCell")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",        &pic::PicCell::dt)
    .def_readwrite("dx",        &pic::PicCell::dx)
    .def_readwrite("container", &pic::PicCell::container);



  py::class_<pic::ParticleBlock>(m, "ParticleBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def("reserve",       &pic::ParticleBlock::reserve)
    .def("resizeEM",      &pic::ParticleBlock::resizeEM)
    .def("add_particle",  &pic::ParticleBlock::add_particle)
    .def("loc",          [](pic::ParticleBlock& s, size_t idim) 
        {
          return s.loc(idim); 
        }, py::return_value_policy::reference)
    .def("vel",          [](pic::ParticleBlock& s, size_t idim) 
        {
          return s.vel(idim); 
        }, py::return_value_policy::reference)
    //temporary binding
    .def("ex",          [](pic::ParticleBlock& s) 
        {
          return s.Epart[0];
        }, py::return_value_policy::reference);
    


    py::class_<pic::Pusher>(m, "Pusher")
      .def(py::init<>())
      .def("solve", &pic::Pusher::solve);

    py::class_<pic::ParticleFieldInterpolator>(m, "ParticleFieldInterpolator")
      .def(py::init<>())
      .def("solve", &pic::ParticleFieldInterpolator::solve);








}
