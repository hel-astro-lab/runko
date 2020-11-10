#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


//--------------------------------------------------
// experimental PIC module
  
#include "../pic/tile.h"
#include "../pic/pushers/pusher.h"
#include "../pic/pushers/boris.h"
#include "../pic/pushers/boris_drag.h"
#include "../pic/pushers/boris_rad.h"
#include "../pic/pushers/boris_grav.h"
#include "../pic/pushers/vay.h"

#include "../pic/interpolators/interpolator.h"
#include "../pic/interpolators/linear.h"

#include "../pic/depositers/depositer.h"
#include "../pic/depositers/zigzag.h"

#include "../pic/communicate.h"

#include "../pic/boundaries/wall.h"
#include "../pic/boundaries/piston.h"

#include "../io/writers/writer.h"
#include "../io/writers/pic.h"
#include "../io/snapshots/test_prtcls.h"
#include "../io/snapshots/pic_moments.h"
#include "../io/tasker.h"





namespace pic {


//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<pic::Tile<D>, 
             fields::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<pic::Tile<D>>
             >(m, 
               pyclass_name.c_str(),
               py::multiple_inheritance()
               )
    .def(py::init<int, int, int>())
    .def_readwrite("cfl",       &pic::Tile<D>::cfl)
    .def("get_container",       &pic::Tile<D>::get_container, 
        py::return_value_policy::reference, py::keep_alive<1,0>())
    .def("set_container",       &pic::Tile<D>::set_container)

    .def("check_outgoing_particles",     &pic::Tile<D>::check_outgoing_particles)
    .def("get_incoming_particles",       &pic::Tile<D>::get_incoming_particles)
    .def("delete_transferred_particles", &pic::Tile<D>::delete_transferred_particles)
    .def("pack_outgoing_particles",      &pic::Tile<D>::pack_outgoing_particles)
    .def("pack_all_particles",           &pic::Tile<D>::pack_all_particles)
    .def("unpack_incoming_particles",    &pic::Tile<D>::unpack_incoming_particles)
    .def("delete_all_particles",         &pic::Tile<D>::delete_all_particles)
    .def("shrink_to_fit_all_particles",  &pic::Tile<D>::shrink_to_fit_all_particles);
}

template<size_t D>
auto declare_prtcl_container(
    py::module& m,
    const std::string& pyclass_name) 
{

  return py::class_<
    pic::ParticleContainer<D>>(m, pyclass_name.c_str())
    .def(py::init<>())
    .def_readwrite("q",   &pic::ParticleContainer<D>::q)
    .def_readwrite("m",   &pic::ParticleContainer<D>::m)
    .def("reserve",       &pic::ParticleContainer<D>::reserve)
    .def("size",          &pic::ParticleContainer<D>::size)
    .def("add_particle",  &pic::ParticleContainer<D>::add_particle)
    .def("add_particle2", [](pic::ParticleContainer<D>& s, 
                            real_prtcl xx, real_prtcl yy, real_prtcl zz,
                            real_prtcl vx, real_prtcl vy, real_prtcl vz, real_prtcl wgt)
        {
          s.add_particle({xx,yy,zz}, {vx,vy,vz}, wgt);
        })
    .def("set_keygen_state", &pic::ParticleContainer<D>::set_keygen_state)
    .def("loc",          [](pic::ParticleContainer<D>& s, size_t idim) 
        {
          return s.loc(idim); 
        }, py::return_value_policy::reference)
    .def("vel",          [](pic::ParticleContainer<D>& s, size_t idim) 
        {
          return s.vel(idim); 
        }, py::return_value_policy::reference)
    .def("wgt",          [](pic::ParticleContainer<D>& s) 
        {
          return s.wgt(); 
        }, py::return_value_policy::reference)
    .def("id",           [](pic::ParticleContainer<D>& s, size_t idim) 
        {
          return s.id(idim); 
        }, py::return_value_policy::reference)

    //temporary binding; only needed for unit tests
    .def("ex",          [](pic::ParticleContainer<D>& s, int i) 
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Epart[i];
        }, py::return_value_policy::reference)
    .def("ey",          [](pic::ParticleContainer<D>& s, int i) 
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Epart[i + 1*nparts];
        }, py::return_value_policy::reference)
    .def("ez",          [](pic::ParticleContainer<D>& s, int i) 
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Epart[i + 2*nparts];
        }, py::return_value_policy::reference)
    .def("bx",          [](pic::ParticleContainer<D>& s, int i)
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Bpart[i];
        }, py::return_value_policy::reference)
    .def("by",          [](pic::ParticleContainer<D>& s, int i) 
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Bpart[i + 1*nparts];
        }, py::return_value_policy::reference)
    .def("bz",          [](pic::ParticleContainer<D>& s, int i) 
        {
          int nparts = s.size();
          assert(i < nparts);
          return s.Bpart[i + 2*nparts];
        }, py::return_value_policy::reference);

}


namespace wall {
  // generator for wall tile
  template<size_t D, int S>
    auto declare_tile(
        py::module& m,
        const std::string& pyclass_name) 
    {
      return
        py::class_<pic::wall::Tile<D, S>,
                   pic::Tile<D>,
                   fields::damping::Tile<D, S>,
                   fields::Tile<D>,
                   corgi::Tile<D>, 
                   std::shared_ptr<pic::wall::Tile<D,S>>
        >(m, 
          pyclass_name.c_str(),
          py::multiple_inheritance()
          )
    .def(py::init<int, int, int>());
    }
}


//--------------------------------------------------

// pybind macros freak out from commas so we hide them using the "using" statement
//using Pusher1D3V = pic::Pusher<1,3>;
//using Pusher2D3V = pic::Pusher<2,3>;
//using Pusher3D3V = pic::Pusher<3,3>;

//template<size_t D>
//using Pusher3V = pic::Pusher<D,3>;

// trampoline class for pic Pusher
//template<size_t D>
//class PyPusher : public Pusher3V<D>
//{
//  void solve(pic::Tile<D>& tile) override {
//  PYBIND11_OVERLOAD_PURE(
//      void,
//      Pusher3V<D>,
//      solve,
//      tile
//      );
//  }
//};


//--------------------------------------------------
// interpolators
  
template<size_t D>
using Interpolator3V = pic::Interpolator<D,3>;

/// trampoline class for pic field Interpolator
template<size_t D>
class PyInterpolator : public Interpolator3V<D>
{
  //using Interpolator3V<D>::Interpolator;

  void solve( pic::Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Interpolator3V<D>,
      solve,
      tile
      );
  }
};

//--------------------------------------------------
// current depositers

template<size_t D>
using Depositer3V = pic::Depositer<D,3>;

/// trampoline class for pic current Depositer
template<size_t D>
class PyDepositer : public Depositer3V<D>
{
  //using Depositer3V<D>::Depositer;

  void solve( pic::Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Depositer3V<D>,
      solve,
      tile
      );
  }
};




//--------------------------------------------------

// python bindings for plasma classes & functions
void bind_pic(py::module& m_sub)
{

  //--------------------------------------------------
  // 1D bindings
  //py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");


  //--------------------------------------------------
  // 2D bindings
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  auto t2 = pic::declare_tile<2>(m_2d, "Tile");
  auto pc2 =pic::declare_prtcl_container<2>(m_2d, "ParticleContainer");

  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");
  auto t3 = pic::declare_tile<3>(m_3d, "Tile");
  auto pc3 =pic::declare_prtcl_container<3>(m_3d, "ParticleContainer");


  //--------------------------------------------------

  // General pusher interface
  //py::class_< pic::Pusher<2,3>, PyPusher<2> > picpusher2d(m_2d, "Pusher");
  py::class_< pic::Pusher<2,3>> picpusher2d(m_2d, "Pusher");
  picpusher2d
    .def(py::init<>())
    .def("solve", py::overload_cast<pic::Tile<2>&>(     &pic::Pusher<2,3>::solve))
    .def("solve", py::overload_cast<pic::Tile<2>&, int>(&pic::Pusher<2,3>::solve));

  // Boris pusher
  py::class_<pic::BorisPusher<2,3>>(m_2d, "BorisPusher", picpusher2d)
    .def(py::init<>());

  // Boris pusher with drag force
  py::class_<pic::BorisPusherDrag<2,3>>(m_2d, "BorisDragPusher", picpusher2d)
    .def_readwrite("drag", &pic::BorisPusherDrag<2,3>::drag)
    .def_readwrite("temp", &pic::BorisPusherDrag<2,3>::temp)
    .def(py::init<>());
    
  // Boris pusher with radiative pressure force
  py::class_<pic::BorisPusherRad<2,3>>(m_2d, "BorisRadPusher", picpusher2d)
    .def_readwrite("drag",      &pic::BorisPusherRad<2,3>::drag)
    .def_readwrite("beam_locx", &pic::BorisPusherRad<2,3>::beam_locx)
    .def(py::init<>());

  // Boris pusher with additional gravity towards cenx 
  py::class_<pic::BorisPusherGrav<2,3>>(m_2d, "BorisGravPusher", picpusher2d)
    .def_readwrite("g0",   &pic::BorisPusherGrav<2,3>::g0)
    .def_readwrite("cenx", &pic::BorisPusherGrav<2,3>::cenx)
    .def(py::init<>());


  // Vay pusher
  py::class_<pic::VayPusher<2,3>>(m_2d, "VayPusher", picpusher2d)
    .def(py::init<>());


  // 3D version
  py::class_< pic::Pusher<3,3>> picpusher3d(m_3d, "Pusher");
  picpusher3d
    .def(py::init<>())
    //.def("solve", &pic::Pusher<3,3>::solve);
    .def("solve", py::overload_cast<pic::Tile<3>&>(     &pic::Pusher<3,3>::solve))
    .def("solve", py::overload_cast<pic::Tile<3>&, int>(&pic::Pusher<3,3>::solve));

  // Boris pusher
  py::class_<pic::BorisPusher<3,3>>(m_3d, "BorisPusher", picpusher3d)
    .def(py::init<>());
    
  // Boris pusher with drag force
  py::class_<pic::BorisPusherDrag<3,3>>(m_3d, "BorisDragPusher", picpusher3d)
    .def_readwrite("drag", &pic::BorisPusherDrag<3,3>::drag)
    .def_readwrite("temp", &pic::BorisPusherDrag<3,3>::temp)
    .def(py::init<>());
    
  // Boris pusher with radiative pressure force
  py::class_<pic::BorisPusherRad<3,3>>(m_3d, "BorisRadPusher", picpusher3d)
    .def_readwrite("drag",      &pic::BorisPusherRad<3,3>::drag)
    .def_readwrite("beam_locx", &pic::BorisPusherRad<3,3>::beam_locx)
    .def(py::init<>());

  // Boris pusher with additional gravity towards cenx 
  py::class_<pic::BorisPusherGrav<3,3>>(m_3d, "BorisGravPusher", picpusher3d)
    .def_readwrite("g0",   &pic::BorisPusherGrav<3,3>::g0)
    .def_readwrite("cenx", &pic::BorisPusherGrav<3,3>::cenx)
    .def(py::init<>());

  // Vay
  py::class_<pic::VayPusher<3,3>>(m_3d, "VayPusher", picpusher3d)
    .def(py::init<>());


  //--------------------------------------------------

  // General interpolator interface
  py::class_< pic::Interpolator<2,3>, PyInterpolator<2> > picinterp2d(m_2d, "Interpolator");
  picinterp2d
    .def(py::init<>())
    .def("solve", &pic::Interpolator<2,3>::solve);

  // Linear pusher
  py::class_<pic::LinearInterpolator<2,3>>(m_2d, "LinearInterpolator", picinterp2d)
    .def(py::init<>());


  // 3D version
  py::class_< pic::Interpolator<3,3>, PyInterpolator<3> > picinterp3d(m_3d, "Interpolator");
  picinterp3d
    .def(py::init<>())
    .def("solve", &pic::Interpolator<3,3>::solve);

  // Linear pusher
  py::class_<pic::LinearInterpolator<3,3>>(m_3d, "LinearInterpolator", picinterp3d)
    .def(py::init<>());

  //--------------------------------------------------
    
  // General current depositer interface
  py::class_< pic::Depositer<2,3>, PyDepositer<2> > picdeposit2d(m_2d, "Depositer");
  picdeposit2d
    .def(py::init<>())
    .def("solve", &pic::Depositer<2,3>::solve);

  // zigzag depositer
  py::class_<pic::ZigZag<2,3>>(m_2d, "ZigZag", picdeposit2d)
    .def(py::init<>());


  // 3D version
  py::class_< pic::Depositer<3,3>, PyDepositer<3> > picdeposit3d(m_3d, "Depositer");
  picdeposit3d
    .def(py::init<>())
    .def("solve", &pic::Depositer<3,3>::solve);

  // zigzag depositer
  py::class_<pic::ZigZag<3,3>>(m_3d, "ZigZag", picdeposit3d)
    .def(py::init<>());




  //--------------------------------------------------
  //2 D piston
  py::class_<pic::Piston<2>>(m_2d, "Piston")
    .def(py::init<>())
    .def_readwrite("walloc",   &pic::Piston<2>::walloc)
    .def_readwrite("gammawall",&pic::Piston<2>::gammawall)
    .def_readwrite("betawall", &pic::Piston<2>::betawall)
    .def("solve",              &pic::Piston<2>::solve)
    .def("field_bc",           &pic::Piston<2>::field_bc);

  //3D piston
  py::class_<pic::Piston<3>>(m_3d, "Piston")
    .def(py::init<>())
    .def_readwrite("walloc",   &pic::Piston<3>::walloc)
    .def_readwrite("gammawall",&pic::Piston<3>::gammawall)
    .def_readwrite("betawall", &pic::Piston<3>::betawall)
    .def("solve",              &pic::Piston<3>::solve)
    .def("field_bc",           &pic::Piston<3>::field_bc);



  //--------------------------------------------------
  // wall

  auto tw1 = pic::wall::declare_tile<2, -1>(m_2d, "Tile_wall_LX");
  auto tw2 = pic::wall::declare_tile<2, +1>(m_2d, "Tile_wall_RX");
  auto tw3 = pic::wall::declare_tile<2, -2>(m_2d, "Tile_wall_LY");
  auto tw4 = pic::wall::declare_tile<2, +2>(m_2d, "Tile_wall_RY");

  //--------------------------------------------------
  // Quick IO 

  py::class_<h5io::TestPrtclWriter<2>>(m_2d, "TestPrtclWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::TestPrtclWriter<2>::write);

  // 3D test particles
  py::class_<h5io::TestPrtclWriter<3>>(m_3d, "TestPrtclWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::TestPrtclWriter<3>::write);

  //--------------------------------------------------
  // physical moments of distribution

  // 2D
  py::class_<h5io::PicMomentsWriter<2>>(m_2d, "PicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::PicMomentsWriter<2>::write);

  // 3D
  py::class_<h5io::PicMomentsWriter<3>>(m_3d, "PicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::PicMomentsWriter<3>::write);



  //--------------------------------------------------
  // Full IO

  // 2D
  m_2d.def("write_particles",  &pic::write_particles<2>);
  m_2d.def("read_particles",   &pic::read_particles<2>);

  // 3D
  m_3d.def("write_particles",  &pic::write_particles<3>);
  m_3d.def("read_particles",   &pic::read_particles<3>);

  //--------------------------------------------------
  // wall
  auto tw13d = pic::wall::declare_tile<3, -1>(m_3d, "Tile_wall_LX");
  auto tw23d = pic::wall::declare_tile<3, +1>(m_3d, "Tile_wall_RX");
  auto tw33d = pic::wall::declare_tile<3, -2>(m_3d, "Tile_wall_LY");
  auto tw43d = pic::wall::declare_tile<3, +2>(m_3d, "Tile_wall_RY");
  auto tw53d = pic::wall::declare_tile<3, -3>(m_3d, "Tile_wall_LZ");
  auto tw63d = pic::wall::declare_tile<3, +3>(m_3d, "Tile_wall_RZ");

}


}
