#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#include "core/ffe/tile.h"

#include "core/ffe/currents/rffe2.h"
#include "core/ffe/currents/rffe4.h"
#include "core/ffe/currents/ffe2.h"
#include "core/ffe/currents/ffe4.h"
#include "core/ffe/slim_grids.h"


namespace ffe {

//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<ffe::Tile<D>, 
             emf::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<ffe::Tile<D>>
             >(m, 
               pyclass_name.c_str(),
               py::multiple_inheritance()
               )
    .def(py::init<int, int, int>())
    .def_readwrite("cfl",       &ffe::Tile<D>::cfl)
    .def("copy_eb",             &ffe::Tile<D>::copy_eb)
    .def("rk3_update",          &ffe::Tile<D>::rk3_update);

}


// python bindings for plasma classes & functions
void bind_ffe(py::module& m_sub)
{

  // skinny version of the Yee lattice with only (e and b meshes)
  py::class_<ffe::SlimGrids>(m_sub, "SlimGrids")
    .def(py::init<int, int, int>())
    .def_readwrite("ex",   &ffe::SlimGrids::ex)
    .def_readwrite("ey",   &ffe::SlimGrids::ey)
    .def_readwrite("ez",   &ffe::SlimGrids::ez)
    .def_readwrite("bx",   &ffe::SlimGrids::bx)
    .def_readwrite("by",   &ffe::SlimGrids::by)
    .def_readwrite("bz",   &ffe::SlimGrids::bz)
    .def("set_grids",        &ffe::SlimGrids::set_grids)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self +  py::self)
    .def(py::self -  py::self)
    .def(py::self *  float())
    .def(py::self /  float());


  m_sub.def("set_step", [](emf::Grids& gs, ffe::SlimGrids sgs)
      -> void 
      {
        gs.ex = sgs.ex;
        gs.ey = sgs.ey;
        gs.ez = sgs.ez;
        
        gs.bx = sgs.bx;
        gs.by = sgs.by;
        gs.bz = sgs.bz;
     }
  );


  //--------------------------------------------------
  // 1D bindings
  //py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");

  //--------------------------------------------------
  // 2D bindings
  //py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  //auto t2 = ffe::declare_tile<2>(m_2d, "Tile");


  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");
  auto t3 = ffe::declare_tile<3>(m_3d, "Tile");


  //--------------------------------------------------
  // 2D Current solver bindings


  //--------------------------------------------------
  // 3D Current solver bindings
  py::class_< ffe::rFFE2<3> > brffe2(m_3d, "rFFE2");
  brffe2
    .def(py::init<int, int, int>())
    .def("comp_rho",     &ffe::rFFE2<3>::comp_rho)
    .def("push_eb",      &ffe::rFFE2<3>::push_eb)
    .def("add_jperp",    &ffe::rFFE2<3>::add_jperp)
    .def("remove_jpar",  &ffe::rFFE2<3>::remove_jpar)
    .def("limit_e",      &ffe::rFFE2<3>::limit_e);

  py::class_< ffe::rFFE4<3> > brffe4(m_3d, "rFFE4");
  brffe4
    .def(py::init<int, int, int>())
    .def("comp_rho",     &ffe::rFFE4<3>::comp_rho)
    .def("push_eb",      &ffe::rFFE4<3>::push_eb)
    .def("add_jperp",    &ffe::rFFE4<3>::add_jperp)
    .def("remove_jpar",  &ffe::rFFE4<3>::remove_jpar)
    .def("limit_e",      &ffe::rFFE4<3>::limit_e);

  py::class_< ffe::FFE2<3> > bffe2(m_3d, "FFE2");
  bffe2
    .def(py::init<int, int, int>())
    .def_readwrite("eta",    &ffe::FFE2<3>::eta)
    .def_readwrite("reltime",&ffe::FFE2<3>::reltime)
    .def("comp_rho",     &ffe::FFE2<3>::comp_rho)
    .def("push_eb",      &ffe::FFE2<3>::push_eb)
    .def("add_jperp",    &ffe::FFE2<3>::add_jperp)
    .def("add_jpar",     &ffe::FFE2<3>::add_jpar)
    .def("limit_e",      &ffe::FFE2<3>::limit_e)
    .def("add_diffusion",&ffe::FFE2<3>::add_diffusion);


  py::class_< ffe::FFE4<3> > bffe4(m_3d, "FFE4");
  bffe4
    .def(py::init<int, int, int>())
    .def_readwrite("eta",    &ffe::FFE4<3>::eta)
    .def_readwrite("reltime",&ffe::FFE4<3>::reltime)
    .def("comp_rho",     &ffe::FFE4<3>::comp_rho)
    .def("push_eb",      &ffe::FFE4<3>::push_eb)
    .def("add_jperp",    &ffe::FFE4<3>::add_jperp)
    .def("add_jpar",     &ffe::FFE4<3>::add_jpar)
    .def("remove_jpar",  &ffe::FFE4<3>::remove_jpar)
    .def("limit_e",      &ffe::FFE4<3>::limit_e)
    .def("add_diffusion",&ffe::FFE4<3>::add_diffusion);



}

} // end of ns ffe
