#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


//--------------------------------------------------
// experimental PIC module
  
#include "core/pic/tile.h"
#include "core/pic/pushers/pusher.h"
#include "core/pic/pushers/boris.h"
#include "core/pic/pushers/boris_drag.h"
#include "core/pic/pushers/boris_rad.h"
#include "core/pic/pushers/boris_grav.h"
#include "core/pic/pushers/vay.h"
#include "core/pic/pushers/higuera_cary.h"
#include "core/pic/pushers/rgca.h"
#include "core/pic/pushers/pulsar.h"
#include "core/pic/pushers/photon.h"

#include "core/pic/interpolators/interpolator.h"
#include "core/pic/interpolators/linear_1st.h"
#include "core/pic/interpolators/quadratic_2nd.h"
#include "core/pic/interpolators/cubic_3rd.h"
#include "core/pic/interpolators/quartic_4th.h"

#include "core/pic/depositers/depositer.h"
#include "core/pic/depositers/zigzag.h"
#include "core/pic/depositers/zigzag_2nd.h"
#include "core/pic/depositers/zigzag_3rd.h"
#include "core/pic/depositers/zigzag_4th.h"
#include "core/pic/depositers/esikerpov_2nd.h"
#include "core/pic/depositers/esikerpov_4th.h"

#include "core/pic/communicate.h"

#include "core/pic/boundaries/wall.h"
#include "core/pic/boundaries/piston.h"
#include "core/pic/boundaries/piston_z.h"
#include "core/pic/boundaries/star_surface_injector.h"

#include "io/writers/writer.h"
#include "io/writers/pic.h"
#include "io/snapshots/test_prtcls.h"
#include "io/snapshots/pic_moments.h"
#include "io/snapshots/master_only_moments.h"
#include "io/tasker.h"





namespace pic {


//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<pic::Tile<D>, 
             emf::Tile<D>,
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
    .def_readwrite("type",&pic::ParticleContainer<D>::type)
    .def("reserve",       &pic::ParticleContainer<D>::reserve)
    .def("size",          &pic::ParticleContainer<D>::size)
    .def("add_particle",  &pic::ParticleContainer<D>::add_particle)
    .def("add_particle2", [](pic::ParticleContainer<D>& s, 
                            float xx, float yy, float zz,
                            float vx, float vy, float vz, float wgt)
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
        }, py::return_value_policy::reference)
    .def("__getitem__", [](pic::ParticleContainer<D>& s, const py::tuple& indx)
      {
        auto ip = indx[0].cast<int>();
        auto v = indx[1].cast<int>();
        int nparts = s.size();

        if (ip >= nparts) throw py::index_error();
        if ( v >= 12)     throw py::index_error();

        if(v == 0) return s.loc(0, ip);
        if(v == 1) return s.loc(1, ip);
        if(v == 2) return s.loc(2, ip);

        if(v == 3) return s.vel(0, ip);
        if(v == 4) return s.vel(1, ip);
        if(v == 5) return s.vel(2, ip);

        if(v == 6) return s.Epart[ip + 0*nparts];
        if(v == 7) return s.Epart[ip + 1*nparts];
        if(v == 8) return s.Epart[ip + 2*nparts];

        if(v == 9) return s.Bpart[ip + 0*nparts];
        if(v ==10) return s.Bpart[ip + 1*nparts];
        if(v ==11) return s.Bpart[ip + 2*nparts];

        // we should not end up here;
        throw py::index_error();
        return float(0.0);
      })
    .def("__setitem__", [](pic::ParticleContainer<D>& s, const py::tuple& indx, float val)
      {
        auto ip = indx[0].cast<int>();
        auto v = indx[1].cast<int>();
        int nparts = s.size();

        if (ip >= nparts) throw py::index_error();
        if ( v >= 12)     throw py::index_error();

        if(v == 0) s.loc(0, ip) = val;
        if(v == 1) s.loc(1, ip) = val;
        if(v == 2) s.loc(2, ip) = val;

        if(v == 3) s.vel(0, ip) = val;
        if(v == 4) s.vel(1, ip) = val;
        if(v == 5) s.vel(2, ip) = val;

        if(v == 6) s.Epart[ip + 0*nparts] = val;
        if(v == 7) s.Epart[ip + 1*nparts] = val;
        if(v == 8) s.Epart[ip + 2*nparts] = val;

        if(v == 9) s.Bpart[ip + 0*nparts] = val;
        if(v ==10) s.Bpart[ip + 1*nparts] = val;
        if(v ==11) s.Bpart[ip + 2*nparts] = val;

        return;
      });

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
                   emf::damping::Tile<D, S>,
                   emf::Tile<D>,
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
  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");
  auto t1 = pic::declare_tile<1>(m_1d, "Tile");
  auto pc1 =pic::declare_prtcl_container<1>(m_1d, "ParticleContainer");

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
  //--------------------------------------------------
  // 1D version
  py::class_< pic::Pusher<1,3>> picpusher1d(m_1d, "Pusher");
  picpusher1d
    .def(py::init<>())
    .def_readwrite("bx_ext",  &pic::Pusher<1,3>::bx_ext)
    .def_readwrite("by_ext",  &pic::Pusher<1,3>::by_ext)
    .def_readwrite("bz_ext",  &pic::Pusher<1,3>::bz_ext)
    .def_readwrite("ex_ext",  &pic::Pusher<1,3>::ex_ext)
    .def_readwrite("ey_ext",  &pic::Pusher<1,3>::ey_ext)
    .def_readwrite("ez_ext",  &pic::Pusher<1,3>::ez_ext)
    //.def("solve", &pic::Pusher<1,3>::solve);
    //.def("solve", py::overload_cast<pic::Tile<1>&>(     &pic::Pusher<1,3>::solve))
    //.def("solve", py::overload_cast<pic::Tile<1>&, int>(&pic::Pusher<1,3>::solve));
    .def("solve", static_cast<void(pic::Pusher<1,3>::*)(pic::Tile<1>&     )>(&pic::Pusher<1,3>::solve))
    .def("solve", static_cast<void(pic::Pusher<1,3>::*)(pic::Tile<1>&, int)>(&pic::Pusher<1,3>::solve));


  // Boris pusher
  py::class_<pic::BorisPusher<1,3>>(m_1d, "BorisPusher", picpusher1d)
    .def(py::init<>());

  // 1D Higuera-Cary pusher
  py::class_<pic::HigueraCaryPusher<1,3>>(m_1d, "HigueraCaryPusher", picpusher1d)
    .def(py::init<>());
    
  // Photon pusher
  py::class_<pic::PhotonPusher<1,3>>(m_1d, "PhotonPusher", picpusher1d)
    .def(py::init<>());

  //--------------------------------------------------
  // 2D version
  //py::class_< pic::Pusher<2,3>, PyPusher<2> > picpusher2d(m_2d, "Pusher");
  py::class_< pic::Pusher<2,3>> picpusher2d(m_2d, "Pusher");
  picpusher2d
    .def(py::init<>())
    .def_readwrite("bx_ext",  &pic::Pusher<2,3>::bx_ext)
    .def_readwrite("by_ext",  &pic::Pusher<2,3>::by_ext)
    .def_readwrite("bz_ext",  &pic::Pusher<2,3>::bz_ext)
    .def_readwrite("ex_ext",  &pic::Pusher<2,3>::ex_ext)
    .def_readwrite("ey_ext",  &pic::Pusher<2,3>::ey_ext)
    .def_readwrite("ez_ext",  &pic::Pusher<2,3>::ez_ext)
    // NOTE: intel compiler has trouble with new overload_cast since its c++14 feature; hence we use the old & nasty c cast
    //.def("solve", py::overload_cast<pic::Tile<2>&>(     &pic::Pusher<2,3>::solve))
    //.def("solve", py::overload_cast<pic::Tile<2>&, int>(&pic::Pusher<2,3>::solve));
    .def("solve", static_cast<void(pic::Pusher<2,3>::*)(pic::Tile<2>&     )>(&pic::Pusher<2,3>::solve))
    .def("solve", static_cast<void(pic::Pusher<2,3>::*)(pic::Tile<2>&, int)>(&pic::Pusher<2,3>::solve));


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

  // Hiuera-Cary pusher
  py::class_<pic::HigueraCaryPusher<2,3>>(m_2d, "HigueraCaryPusher", picpusher2d)
    .def(py::init<>());

  // reduced guiding center approximation
  py::class_<pic::rGCAPusher<2,3>>(m_2d, "rGCAPusher", picpusher2d)
    .def(py::init<>());

  // special rGCA + gravity pusher for pulsars
  py::class_<pic::PulsarPusher<2,3>>(m_2d, "PulsarPusher", picpusher2d)
    .def_readwrite("radius",    &pic::PulsarPusher<2,3>::rad_star)
    .def_readwrite("radius_pc", &pic::PulsarPusher<2,3>::rad_pcap)
    .def_readwrite("period",    &pic::PulsarPusher<2,3>::period_star)
    .def_readwrite("cenx",      &pic::PulsarPusher<2,3>::cenx)
    .def_readwrite("ceny",      &pic::PulsarPusher<2,3>::ceny)
    .def_readwrite("cenz",      &pic::PulsarPusher<2,3>::cenz)
    .def_readwrite("grav_const",&pic::PulsarPusher<2,3>::gravity_const)
    //.def_readwrite("B0",       &pic::PulsarPusher<2>::B0)
    //.def_readwrite("chi",      &pic::PulsarPusher<2>::chi)
    //.def_readwrite("phase",    &pic::PulsarPusher<2>::phase)
    //.def_readwrite("delta",    &pic::PulsarPusher<2>::delta)
    //.def_readwrite("delta_pc", &pic::PulsarPusher<2>::delta_pc)
    .def(py::init<>());

  // Photon pusher
  py::class_<pic::PhotonPusher<2,3>>(m_2d, "PhotonPusher", picpusher2d)
    .def(py::init<>());

   //--------------------------------------------------
  // 3D version
  py::class_< pic::Pusher<3,3>> picpusher3d(m_3d, "Pusher");
  picpusher3d
    .def(py::init<>())
    .def_readwrite("bx_ext",  &pic::Pusher<3,3>::bx_ext)
    .def_readwrite("by_ext",  &pic::Pusher<3,3>::by_ext)
    .def_readwrite("bz_ext",  &pic::Pusher<3,3>::bz_ext)
    .def_readwrite("ex_ext",  &pic::Pusher<3,3>::ex_ext)
    .def_readwrite("ey_ext",  &pic::Pusher<3,3>::ey_ext)
    .def_readwrite("ez_ext",  &pic::Pusher<3,3>::ez_ext)
    //.def("solve", &pic::Pusher<3,3>::solve);
    //.def("solve", py::overload_cast<pic::Tile<3>&>(     &pic::Pusher<3,3>::solve))
    //.def("solve", py::overload_cast<pic::Tile<3>&, int>(&pic::Pusher<3,3>::solve));
    .def("solve", static_cast<void(pic::Pusher<3,3>::*)(pic::Tile<3>&     )>(&pic::Pusher<3,3>::solve))
    .def("solve", static_cast<void(pic::Pusher<3,3>::*)(pic::Tile<3>&, int)>(&pic::Pusher<3,3>::solve));


  // Boris pusher
  py::class_<pic::BorisPusher<3,3>>(m_3d, "BorisPusher", picpusher3d)
    .def(py::init<>());
    
  // Boris pusher with drag force
  py::class_<pic::BorisPusherDrag<3,3>>(m_3d, "BorisDragPusher", picpusher3d)
    .def_readwrite("drag", &pic::BorisPusherDrag<3,3>::drag)
    .def_readwrite("temp", &pic::BorisPusherDrag<3,3>::temp)
    .def_readwrite("freezing_factor", &pic::BorisPusherDrag<3,3>::freezing_factor)
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

  // Higuera-Cary
  py::class_<pic::HigueraCaryPusher<3,3>>(m_3d, "HigueraCaryPusher", picpusher3d)
    .def(py::init<>());

  // reduced guiding center approximation
  py::class_<pic::rGCAPusher<3,3>>(m_3d, "rGCAPusher", picpusher3d)
    .def(py::init<>());

  // special rGCA + gravity pusher for pulsars
  py::class_<pic::PulsarPusher<3,3>>(m_3d, "PulsarPusher", picpusher3d)
    .def_readwrite("radius",    &pic::PulsarPusher<3,3>::rad_star)
    .def_readwrite("radius_pc", &pic::PulsarPusher<3,3>::rad_pcap)
    .def_readwrite("period",    &pic::PulsarPusher<3,3>::period_star)
    .def_readwrite("cenx",      &pic::PulsarPusher<3,3>::cenx)
    .def_readwrite("ceny",      &pic::PulsarPusher<3,3>::ceny)
    .def_readwrite("cenz",      &pic::PulsarPusher<3,3>::cenz)
    .def_readwrite("grav_const",&pic::PulsarPusher<3,3>::gravity_const)
    //.def_readwrite("B0",       &pic::PulsarPusher<3>::B0)
    //.def_readwrite("chi",      &pic::PulsarPusher<3>::chi)
    //.def_readwrite("phase",    &pic::PulsarPusher<3>::phase)
    //.def_readwrite("delta",    &pic::PulsarPusher<3>::delta)
    //.def_readwrite("delta_pc", &pic::PulsarPusher<3>::delta_pc)
    .def(py::init<>());

  // Photon pusher
  py::class_<pic::PhotonPusher<3,3>>(m_3d, "PhotonPusher", picpusher3d)
    .def(py::init<>());

  //--------------------------------------------------

  // General interpolator interface
  //--------------------------------------------------
  // 1D version
  py::class_< pic::Interpolator<1,3>, PyInterpolator<1> > picinterp1d(m_1d, "Interpolator");
  picinterp1d
    .def(py::init<>())
    .def("solve", &pic::Interpolator<1,3>::solve);

  // Linear pusher
  py::class_<pic::LinearInterpolator<1,3>>(m_1d, "LinearInterpolator", picinterp1d)
    .def(py::init<>());

  //--------------------------------------------------
  // 2D version
  py::class_< pic::Interpolator<2,3>, PyInterpolator<2> > picinterp2d(m_2d, "Interpolator");
  picinterp2d
    .def(py::init<>())
    .def("solve", &pic::Interpolator<2,3>::solve);

  // Linear pusher
  py::class_<pic::LinearInterpolator<2,3>>(m_2d, "LinearInterpolator", picinterp2d)
    .def(py::init<>());

  // 2nd order quadratic
  py::class_<pic::QuadraticInterpolator<2>>(m_2d, "QuadraticInterpolator", picinterp2d)
    .def(py::init<>());


  //--------------------------------------------------
  // 3D version
  py::class_< pic::Interpolator<3,3>, PyInterpolator<3> > picinterp3d(m_3d, "Interpolator");
  picinterp3d
    .def(py::init<>())
    .def("solve", &pic::Interpolator<3,3>::solve);

  // Linear pusher
  py::class_<pic::LinearInterpolator<3,3>>(m_3d, "LinearInterpolator", picinterp3d)
    .def(py::init<>());

  // 2nd order quadratic
  py::class_<pic::QuadraticInterpolator<3>>(m_3d, "QuadraticInterpolator", picinterp3d)
    .def(py::init<>());

  // 3rd order cubic
  py::class_<pic::CubicInterpolator<3>>(m_3d, "CubicInterpolator", picinterp3d)
    .def(py::init<>());

  // 4th order quartic
  py::class_<pic::QuarticInterpolator<3>>(m_3d, "QuarticInterpolator", picinterp3d)
    .def(py::init<>());

  //--------------------------------------------------
    
  // General current depositer interface
  //--------------------------------------------------
  // 1D version
  py::class_< pic::Depositer<1,3>, PyDepositer<1> > picdeposit1d(m_1d, "Depositer");
  picdeposit1d
    .def(py::init<>())
    .def("solve", &pic::Depositer<1,3>::solve);

  // zigzag depositer
  py::class_<pic::ZigZag<1,3>>(m_1d, "ZigZag", picdeposit1d)
    .def(py::init<>());

  //--------------------------------------------------
  // 2D version
  py::class_< pic::Depositer<2,3>, PyDepositer<2> > picdeposit2d(m_2d, "Depositer");
  picdeposit2d
    .def(py::init<>())
    .def("solve", &pic::Depositer<2,3>::solve);

  // zigzag depositer
  py::class_<pic::ZigZag<2,3>>(m_2d, "ZigZag", picdeposit2d)
    .def(py::init<>());

  py::class_<pic::ZigZag_2nd<2,3>>(m_2d, "ZigZag_2nd", picdeposit2d)
    .def(py::init<>());

  py::class_<pic::ZigZag_3rd<2,3>>(m_2d, "ZigZag_3rd", picdeposit2d)
    .def(py::init<>());

  py::class_<pic::ZigZag_4th<2,3>>(m_2d, "ZigZag_4th", picdeposit2d)
    .def(py::init<>());


  //--------------------------------------------------
  // 3D version
  py::class_< pic::Depositer<3,3>, PyDepositer<3> > picdeposit3d(m_3d, "Depositer");
  picdeposit3d
    .def(py::init<>())
    .def("solve", &pic::Depositer<3,3>::solve);

  // zigzag depositer
  py::class_<pic::ZigZag<3,3>>(m_3d, "ZigZag", picdeposit3d)
    .def(py::init<>());

  py::class_<pic::ZigZag_2nd<3,3>>(m_3d, "ZigZag_2nd", picdeposit3d)
    .def(py::init<>());

  py::class_<pic::ZigZag_3rd<3,3>>(m_3d, "ZigZag_3rd", picdeposit3d)
    .def(py::init<>());

  py::class_<pic::ZigZag_4th<3,3>>(m_3d, "ZigZag_4th", picdeposit3d)
    .def(py::init<>());

  py::class_<pic::Esikerpov_2nd<3,3>>(m_3d, "Esikerpov_2nd", picdeposit3d)
    .def(py::init<>());

  py::class_<pic::Esikerpov_4th<3,3>>(m_3d, "Esikerpov_4th", picdeposit3d)
    .def(py::init<>());

  //--------------------------------------------------

  //1 D piston
  py::class_<pic::Piston<1>>(m_1d, "Piston")
    .def(py::init<>())
    .def_readwrite("walloc",   &pic::Piston<1>::walloc)
    .def_readwrite("gammawall",&pic::Piston<1>::gammawall)
    .def_readwrite("betawall", &pic::Piston<1>::betawall)
    .def("solve",              &pic::Piston<1>::solve)
    .def("field_bc",           &pic::Piston<1>::field_bc);

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

  //3D piston
  py::class_<pic::PistonZdir<3>>(m_3d, "PistonZdir")
    .def(py::init<>())
    .def_readwrite("walloc",   &pic::PistonZdir<3>::wallocz)
    .def_readwrite("wdir",     &pic::PistonZdir<3>::wdir)
    .def("solve",              &pic::PistonZdir<3>::solve)
    .def("clean_prtcls",       &pic::PistonZdir<3>::clean_prtcls)
    .def("field_bc",           &pic::PistonZdir<3>::field_bc);


  //--------------------------------------------------
  // 1D wall

  auto tw11d = pic::wall::declare_tile<1, -1>(m_1d, "Tile_wall_LX");
  auto tw21d = pic::wall::declare_tile<1, +1>(m_1d, "Tile_wall_RX");

  //--------------------------------------------------
  // 2D wall

  auto tw12d = pic::wall::declare_tile<2, -1>(m_2d, "Tile_wall_LX");
  auto tw22d = pic::wall::declare_tile<2, +1>(m_2d, "Tile_wall_RX");
  auto tw32d = pic::wall::declare_tile<2, -2>(m_2d, "Tile_wall_LY");
  auto tw42d = pic::wall::declare_tile<2, +2>(m_2d, "Tile_wall_RY");

  //--------------------------------------------------
  // Quick IO 

  // 1D test particles
  py::class_<h5io::TestPrtclWriter<1>>(m_1d, "TestPrtclWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::TestPrtclWriter<1>::write)
    .def_readwrite("ispc", &h5io::TestPrtclWriter<1>::ispc);
  
  // 2D test particles
  py::class_<h5io::TestPrtclWriter<2>>(m_2d, "TestPrtclWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::TestPrtclWriter<2>::write)
    .def_readwrite("ispc", &h5io::TestPrtclWriter<2>::ispc);

  // 3D test particles
  py::class_<h5io::TestPrtclWriter<3>>(m_3d, "TestPrtclWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::TestPrtclWriter<3>::write)
    .def_readwrite("ispc", &h5io::TestPrtclWriter<3>::ispc);

  //--------------------------------------------------
  // physical moments of distribution

  // 1D
  py::class_<h5io::PicMomentsWriter<1>>(m_1d, "PicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::PicMomentsWriter<1>::write);
  
  // 2D
  py::class_<h5io::PicMomentsWriter<2>>(m_2d, "PicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",       &h5io::PicMomentsWriter<2>::write)
    .def("get_slice", [](h5io::PicMomentsWriter<2> &s, int k)
            {
                const auto nx = static_cast<pybind11::ssize_t>( s.nx );
                const auto ny = static_cast<pybind11::ssize_t>( s.ny );
                auto v = pybind11::array_t<float>( {nx, ny}, s.arrs[k].data() );
                return v;
            });

  // 3D
  py::class_<h5io::PicMomentsWriter<3>>(m_3d, "PicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::PicMomentsWriter<3>::write);

  // 3D
  py::class_<h5io::MasterPicMomentsWriter<3>>(m_3d, "MasterPicMomentsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write", &h5io::MasterPicMomentsWriter<3>::write);


  //--------------------------------------------------
  // Full IO

  // 1D
  m_1d.def("write_particles",  &pic::write_particles<1>);
  m_1d.def("read_particles",   &pic::read_particles<1>);
  
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

  //--------------------------------------------------
  // star


  // 2D rotating conductor
  py::class_<pic::Star<2>>(m_2d, "Star")
    .def(py::init<>())
    .def_readwrite("B0",             &pic::Star<2>::B0)
    .def_readwrite("radius",         &pic::Star<2>::radius)
    .def_readwrite("period",         &pic::Star<2>::period)
    .def_readwrite("chi_mu",         &pic::Star<2>::chi_mu)
    .def_readwrite("chi_om",         &pic::Star<2>::chi_om)
    .def_readwrite("phase_mu",       &pic::Star<2>::phase_mu)
    .def_readwrite("phase_om",       &pic::Star<2>::phase_om)
    .def_readwrite("cenx",           &pic::Star<2>::cenx)
    .def_readwrite("ceny",           &pic::Star<2>::ceny)
    .def_readwrite("cenz",           &pic::Star<2>::cenz)
    .def_readwrite("delta",          &pic::Star<2>::delta)
    .def_readwrite("radius_pc",      &pic::Star<2>::radius_pc)
    .def_readwrite("delta_pc",       &pic::Star<2>::delta_pc)
    .def_readwrite("Nx",             &pic::Star<2>::Nx)
    .def_readwrite("Ny",             &pic::Star<2>::Ny)
    .def_readwrite("Nz",             &pic::Star<2>::Nz)
    .def_readwrite("temp_pairs",     &pic::Star<2>::temp_pairs)
    .def_readwrite("temp_phots",     &pic::Star<2>::temp_phots)
    .def_readwrite("ninj_pairs",     &pic::Star<2>::ninj_pairs)
    .def_readwrite("ninj_phots",     &pic::Star<2>::ninj_phots)
    .def_readwrite("ninj_min_pairs", &pic::Star<2>::ninj_min_pairs)
    .def_readwrite("ninj_min_phots", &pic::Star<2>::ninj_min_phots)
    .def("insert_em",                &pic::Star<2>::insert_em)
    .def("update_b",                 &pic::Star<2>::update_b)
    .def("update_e",                 &pic::Star<2>::update_e)
    .def("solve",                    &pic::Star<2>::solve);


  // 3D rotating conductor
  py::class_<pic::Star<3>>(m_3d, "Star")
    .def(py::init<>())
    .def_readwrite("B0",             &pic::Star<3>::B0)
    .def_readwrite("radius",         &pic::Star<3>::radius)
    .def_readwrite("period",         &pic::Star<3>::period)
    .def_readwrite("chi_mu",         &pic::Star<3>::chi_mu)
    .def_readwrite("chi_om",         &pic::Star<3>::chi_om)
    .def_readwrite("phase_mu",       &pic::Star<3>::phase_mu)
    .def_readwrite("phase_om",       &pic::Star<3>::phase_om)
    .def_readwrite("cenx",           &pic::Star<3>::cenx)
    .def_readwrite("ceny",           &pic::Star<3>::ceny)
    .def_readwrite("cenz",           &pic::Star<3>::cenz)
    .def_readwrite("delta",          &pic::Star<3>::delta)
    .def_readwrite("radius_pc",      &pic::Star<3>::radius_pc)
    .def_readwrite("delta_pc",       &pic::Star<3>::delta_pc)
    .def_readwrite("Nx",             &pic::Star<3>::Nx)
    .def_readwrite("Ny",             &pic::Star<3>::Ny)
    .def_readwrite("Nz",             &pic::Star<3>::Nz)
    .def_readwrite("temp_pairs",     &pic::Star<3>::temp_pairs)
    .def_readwrite("temp_phots",     &pic::Star<3>::temp_phots)
    .def_readwrite("ninj_pairs",     &pic::Star<3>::ninj_pairs)
    .def_readwrite("ninj_phots",     &pic::Star<3>::ninj_phots)
    .def_readwrite("ninj_min_pairs", &pic::Star<3>::ninj_min_pairs)
    .def_readwrite("ninj_min_phots", &pic::Star<3>::ninj_min_phots)
    .def("insert_em",                &pic::Star<3>::insert_em)
    .def("update_b",                 &pic::Star<3>::update_b)
    .def("update_e",                 &pic::Star<3>::update_e)
    .def("solve",                    &pic::Star<3>::solve);


}


}
