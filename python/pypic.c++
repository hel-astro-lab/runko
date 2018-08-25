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

#include "../pic/boundaries/wall.h"


namespace pic {


//--------------------------------------------------

namespace wall {
  // generator for wall tile
  template<size_t D, int S>
    auto declare_Tile(
        py::module& m,
        const std::string& pyclass_name) 
    {
      return
        py::class_<pic::wall::Tile<D, S>,
      fields::damping::Tile<D, S>,
      pic::Tile<D>,
      std::shared_ptr<pic::wall::Tile<D,S>>
        >(m, pyclass_name.c_str() );
    }
}



// python bindings for plasma classes & functions
void bind_pic(py::module& m)
{

  // Loading tile bindings from corgi library
  //py::object corgiTile<2> = (py::object) py::module::import("pycorgi").attr("Tile<2>");
  //py::object plasmaTile<2> = (py::object) py::module::import("pyplasma").attr("Tile<2>");


  py::class_<pic::Tile<2>, 
             fields::Tile<2>,
             corgi::Tile<2>, 
             std::shared_ptr<pic::Tile<2>>
             >(m, "Tile2D")
    .def(py::init<size_t, size_t>())
    .def_readwrite("dx",        &pic::Tile<2>::dx)
    .def_readwrite("cfl",       &pic::Tile<2>::cfl)
    //.def_readwrite("container", &pic::Tile<2>::container);
    .def("get_container",       &pic::Tile<2>::get_container, 
        py::return_value_policy::reference)
    .def("set_container",       &pic::Tile<2>::set_container);



  py::class_<pic::ParticleBlock>(m, "ParticleBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("q",   &pic::ParticleBlock::q)
    .def("reserve",       &pic::ParticleBlock::reserve)
    .def("resizeEM",      &pic::ParticleBlock::resizeEM)
    .def("add_particle",  &pic::ParticleBlock::add_particle)
    .def("add_particle2", [](pic::ParticleBlock& s, 
                            Realf xx, Realf yy, Realf zz,
                            Realf vx, Realf vy, Realf vz)
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

  //--------------------------------------------------
  // wall

  auto tw1 = pic::wall::declare_Tile<2, -1>(m, "Tile2D_wall_LX");
  auto tw2 = pic::wall::declare_Tile<2, +1>(m, "Tile2D_wall_RX");
  auto tw3 = pic::wall::declare_Tile<2, -2>(m, "Tile2D_wall_LY");
  auto tw4 = pic::wall::declare_Tile<2, +2>(m, "Tile2D_wall_RY");

  tw1.def(py::init<size_t, size_t>());
  tw2.def(py::init<size_t, size_t>());
  tw3.def(py::init<size_t, size_t>());
  tw4.def(py::init<size_t, size_t>());



}


}
