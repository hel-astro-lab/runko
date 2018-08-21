#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


//--------------------------------------------------
// experimental PIC module
  
#include "../pic/tile.h"
#include "../pic/pusher.h"
#include "../pic/field_interpolator.h"
#include "../pic/communicate.h"
#include "../pic/current_deposit.h"
#include "../pic/analyzer.h"
#include "../pic/filters.h"




//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(pypic, m) {

  // Loading tile bindings from corgi library
  py::object corgiTile = (py::object) py::module::import("pycorgi").attr("Tile");

  py::object plasmaTile = (py::object) py::module::import("pyplasma").attr("PlasmaTile");


  py::class_<pic::PicTile, 
             fields::PlasmaTile,
             corgi::Tile, 
             std::shared_ptr<pic::PicTile>
             >(m, "PicTile")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",        &pic::PicTile::dt)
    .def_readwrite("dx",        &pic::PicTile::dx)
    .def_readwrite("cfl",       &pic::PicTile::cfl)
    //.def_readwrite("container", &pic::PicTile::container);
    .def("get_container",       &pic::PicTile::get_container, 
        py::return_value_policy::reference)
    .def("set_container",       &pic::PicTile::set_container);



  py::class_<pic::ParticleBlock>(m, "ParticleBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("q",  &pic::ParticleBlock::q)
    .def("reserve",       &pic::ParticleBlock::reserve)
    .def("resizeEM",      &pic::ParticleBlock::resizeEM)
    .def("add_particle",  &pic::ParticleBlock::add_particle)
    .def("add_particle2", [](pic::ParticleBlock& s, 
                            double xx, double yy, double zz,
                            double vx, double vy, double vz)
        {
          s.add_particle({xx,yy,zz}, {vx,vy,vz});
        })
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
        }, py::return_value_policy::reference)
    .def("ey",          [](pic::ParticleBlock& s) 
        {
          return s.Epart[1];
        }, py::return_value_policy::reference)
    .def("ez",          [](pic::ParticleBlock& s) 
        {
          return s.Epart[2];
        }, py::return_value_policy::reference)
    .def("bx",          [](pic::ParticleBlock& s) 
        {
          return s.Bpart[0];
        }, py::return_value_policy::reference)
    .def("by",          [](pic::ParticleBlock& s) 
        {
          return s.Bpart[1];
        }, py::return_value_policy::reference)
    .def("bz",          [](pic::ParticleBlock& s) 
        {
          return s.Bpart[2];
        }, py::return_value_policy::reference);
    


    py::class_<pic::Pusher>(m, "Pusher")
      .def(py::init<>())
      .def("solve", &pic::Pusher::solve);

    py::class_<pic::ParticleFieldInterpolator>(m, "ParticleFieldInterpolator")
      .def(py::init<>())
      .def("solve", &pic::ParticleFieldInterpolator::solve);


    py::class_<pic::Communicator>(m, "Communicator")
      .def(py::init<>())
      .def("check_outgoing_particles",    &pic::Communicator::check_outgoing_particles)
      .def("get_incoming_particles",      &pic::Communicator::get_incoming_particles)
      .def("delete_transferred_particles",&pic::Communicator::delete_transferred_particles);


    py::class_<pic::Depositer>(m, "Depositer")
      .def(py::init<>())
      .def("deposit", &pic::Depositer::deposit);


    /// Pic tile analyzator
    py::class_<pic::Analyzator>(m, "Analyzator")
      .def(py::init<>())
      .def("analyze", &pic::Analyzator::analyze);


    py::class_<pic::Filter>(m, "Filter")
      .def(py::init<int, int>())
      .def("init_kernel",         &pic::Filter::init_kernel)
      .def("init_gaussian_kernel",&pic::Filter::init_gaussian_kernel)
      //.def("init_sinc_kernel",    &pic::Filter::init_sinc_kernel)
      .def("init_lowpass_fft_kernel",&pic::Filter::init_lowpass_fft_kernel)
      .def("init_3point",            &pic::Filter::init_3point_kernel)
      .def("fft_kernel",             &pic::Filter::fft_kernel)
      .def("fft_image_forward",      &pic::Filter::fft_image_forward)
      .def("fft_image_backward",     &pic::Filter::fft_image_backward)
      .def("apply_kernel",           &pic::Filter::apply_kernel)
      .def("get_padded_current",     &pic::Filter::get_padded_current)
      .def("set_current",            &pic::Filter::set_current)
      .def("direct_convolve_3point", &pic::Filter::direct_convolve_3point)
      .def("set_image",              &pic::Filter::set_image)
      .def("set_kernel",             &pic::Filter::set_kernel)
      .def("get_kernel",             &pic::Filter::get_kernel, py::return_value_policy::reference)
      .def("get_image",              &pic::Filter::get_image,  py::return_value_policy::reference);



}
