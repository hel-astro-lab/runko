#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


//--------------------------------------------------
// experimental PIC module
  
#include "../pic/tile.h"
#include "../pic/solvers/pusher.h"
#include "../pic/solvers/boris.h"

#include "../pic/field_interpolator.h"
#include "../pic/communicate.h"
#include "../pic/current_deposit.h"
#include "../pic/analyzer.h"
#include "../pic/filters.h"

#include "../pic/boundaries/wall.h"


namespace pic {


//--------------------------------------------------
template<size_t D>
auto declare_Tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<pic::Tile<D>, 
             fields::Tile<D>,
             std::shared_ptr<pic::Tile<D>>
             >(m, pyclass_name.c_str())
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("dx",        &pic::Tile<D>::dx)
    .def_readwrite("cfl",       &pic::Tile<D>::cfl)
    //.def_readwrite("container", &pic::Tile<D>::container);
    .def("get_container",       &pic::Tile<D>::get_container, 
        py::return_value_policy::reference)
    .def("set_container",       &pic::Tile<D>::set_container);
}


namespace wall {
  // generator for wall tile
  template<size_t D, int S>
    auto declare_Tile(
        py::module& m,
        const std::string& pyclass_name) 
    {
      return
        py::class_<pic::wall::Tile<D, S>,
      pic::Tile<D>,
      fields::damping::Tile<D, S>,
      std::shared_ptr<pic::wall::Tile<D,S>>
        >(m, pyclass_name.c_str() )
    .def(py::init<size_t, size_t, size_t>());

    }
}


// pybind macros freak out from commas so we hide them using the "using" statement
//using Pusher1D3V = pic::Pusher<1,3>;
//using Pusher2D3V = pic::Pusher<2,3>;
//using Pusher3D3V = pic::Pusher<3,3>;

template<size_t D>
using Pusher3V = pic::Pusher<D,3>;

/// trampoline class for pic Pusher
template<size_t D>
class PyPusher : public Pusher3V<D>
{
  //using Pusher3V<D>::Pusher;

  void solve( pic::Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Pusher3V<D>,
      solve,
      tile
      );
  }
};





// python bindings for plasma classes & functions
void bind_pic(py::module& m_sub)
{


  py::class_<pic::ParticleBlock>(m_sub, "ParticleBlock")
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
    

  //--------------------------------------------------
  // 1D bindings
  //py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");


  //--------------------------------------------------
  // 2D bindings
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");

  auto t2 = pic::declare_Tile<2>(m_2d, "Tile");


  //--------------------------------------------------
  //py::class_<pic::Pusher>(m_2d, "Pusher")
  //  .def(py::init<>())
  //  .def("solve", &pic::Pusher::solve);

  // General pusher interface
  py::class_< pic::Pusher<2,3>, PyPusher<2> > picpusher2d(m_2d, "Pusher");
  picpusher2d
    .def(py::init<>())
    .def("solve", &pic::Pusher<2,3>::solve);

  // Boris pusher
  py::class_<pic::BorisPusher<2,3>>(m_2d, "BorisPusher", picpusher2d)
    .def(py::init<>());





  py::class_<pic::ParticleFieldInterpolator>(m_2d, "ParticleFieldInterpolator")
    .def(py::init<>())
    .def("solve", &pic::ParticleFieldInterpolator::solve);


    py::class_<pic::Communicator>(m_2d, "Communicator")
      .def(py::init<>())
      .def("check_outgoing_particles",    &pic::Communicator::check_outgoing_particles)
      .def("get_incoming_particles",      &pic::Communicator::get_incoming_particles)
      .def("delete_transferred_particles",&pic::Communicator::delete_transferred_particles);


    py::class_<pic::Depositer>(m_2d, "Depositer")
      .def(py::init<>())
      .def("deposit", &pic::Depositer::deposit);


    /// Pic tile analyzator
    py::class_<pic::Analyzator>(m_2d, "Analyzator")
      .def(py::init<>())
      .def("analyze", &pic::Analyzator::analyze);


    py::class_<pic::Filter>(m_2d, "Filter")
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

  auto tw1 = pic::wall::declare_Tile<2, -1>(m_2d, "Tile_wall_LX");
  auto tw2 = pic::wall::declare_Tile<2, +1>(m_2d, "Tile_wall_RX");
  auto tw3 = pic::wall::declare_Tile<2, -2>(m_2d, "Tile_wall_LY");
  auto tw4 = pic::wall::declare_Tile<2, +2>(m_2d, "Tile_wall_RY");


}


}
