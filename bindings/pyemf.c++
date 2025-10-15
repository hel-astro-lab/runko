#include <string>
#include <pybind11/numpy.h>

#include "py_submodules.h"

#include "external/corgi/internals.h"

#include "definitions.h"
#include "tools/mesh.h"

#include "core/emf/tile.h"

#include "core/emf/propagators/propagator.h"
#include "core/emf/propagators/fdtd2.h"
#include "core/emf/propagators/fdtd2_pml.h"
#include "core/emf/propagators/fdtd4.h"
#include "core/emf/propagators/fdtd_general.h"

#include "core/emf/filters/filter.h"
#include "core/emf/filters/binomial2.h"
#include "core/emf/filters/compensator.h"
#include "core/emf/filters/strided_binomial.h"
#include "core/emf/filters/general_binomial.h"
#include "core/emf/filters/sweeping_binomial.h"


#include "core/emf/boundaries/damping_tile.h"
#include "core/emf/boundaries/conductor.h"

#include "io/writers/writer.h"
#include "io/writers/fields.h"
#include "io/snapshots/fields.h"
#include "io/snapshots/master_only_fields.h"
#include "io/snapshots/field_slices.h"
#include "io/tasker.h"


//--------------------------------------------------
  
namespace emf{

  namespace py = pybind11;

  
// generator for damped tile for various directions
template<int D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{
  std::vector<int> iarr{0,1,2};

  return py::class_<
              emf::Tile<D>,
              corgi::Tile<D>, 
              std::shared_ptr<emf::Tile<D>>
            >(m, pyclass_name.c_str() )
    .def(py::init<int, int, int>())
    .def_readwrite("cfl",       &emf::Tile<D>::cfl)
    .def("clear_current",       &emf::Tile<D>::clear_current)
    .def("deposit_current",     &emf::Tile<D>::deposit_current)
    .def("exchange_currents",   &emf::Tile<D>::exchange_currents)
    .def("update_boundaries",   &emf::Tile<D>::update_boundaries,
            py::arg("grid"),
            py::arg("iarr")=iarr)
    .def("get_grids",             &emf::Tile<D>::get_grids,
        py::arg("i")=0,
        py::return_value_policy::reference,
        // keep alive for the lifetime of the grid
        //
        // pybind11:
        // argument indices start at one, while zero refers to the return 
        // value. For methods, index one refers to the implicit this pointer, 
        // while regular arguments begin at index two. 
        // py::keep_alive<nurse,patient>()
        py::keep_alive<1,0>()
        );
}



// generator for damped tile for various directions
template<int D, int S>
auto declare_TileDamped(
    py::module& m,
    const std::string& pyclass_name) 
{
  // using Class = emf::TileDamped<D,S>; 
  // does not function properly; maybe not triggering template?
  // have to use explicit name instead like this

  return py::class_<
             emf::damping::Tile<D,S>,
             emf::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<emf::damping::Tile<D,S>>
          >(m, 
            pyclass_name.c_str(),
            py::multiple_inheritance()
            )
  .def(py::init<int, int, int>())
  .def_readwrite("ex_ref",   &emf::damping::Tile<D,S>::ex_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("ey_ref",   &emf::damping::Tile<D,S>::ey_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("ez_ref",   &emf::damping::Tile<D,S>::ez_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("bx_ref",   &emf::damping::Tile<D,S>::bx_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("by_ref",   &emf::damping::Tile<D,S>::by_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("bz_ref",   &emf::damping::Tile<D,S>::bz_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("fld1",     &emf::damping::Tile<D,S>::fld1)
  .def_readwrite("fld2",     &emf::damping::Tile<D,S>::fld2)
  .def_readwrite("ksupp",    &emf::damping::Tile<D,S>::ksupp)
  .def("damp_fields",        &emf::damping::Tile<D,S>::damp_fields);
}


/// trampoline class for emf Propagator
template<int D>
class PyPropagator : public Propagator<D>
{
  using Propagator<D>::Propagator;

  void push_e( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Propagator<D>,
      push_e,
      tile
      );
  }

  void push_half_b( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Propagator<D>,
      push_half_b,
      tile
      );
  }
};


template<size_t D>
class PyFDTD2 : public FDTD2<D>
{
  using FDTD2<D>::FDTD2;

  void push_e( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTD2<D>, push_e, tile);
  }

  void push_half_b( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTD2<D>, push_half_b, tile);
  }
};

template<size_t D>
class PyFDTD4 : public FDTD4<D>
{
  using FDTD4<D>::FDTD4;

  void push_e( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTD4<D>, push_e, tile);
  }

  void push_half_b( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTD4<D>, push_half_b, tile);
  }
};


template<size_t D>
class PyFDTDGen : public FDTDGen<D>
{
  using FDTDGen<D>::FDTDGen;

  void push_e( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTDGen<D>, push_e, tile);
  }

  void push_half_b( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE( void, FDTDGen<D>, push_half_b, tile);
  }
};

// templated base-class for Propagator; 
// see https://github.com/pybind/pybind11/blob/master/tests/test_virtual_functions.cpp
//template <class Base = Propagator<1>, size_t D=1>
//class PyPropagator : public Base
//{
//  using Base::Base; // inherit constructor
//
//  void push_e( Tile<D>& tile ) override {
//  PYBIND11_OVERLOAD_PURE(
//      void,
//      Base,
//      push_e,
//      tile
//      );
//  }
//
//  void push_half_b( Tile<D>& tile ) override {
//  PYBIND11_OVERLOAD_PURE(
//      void,
//      Base,
//      push_half_b,
//      tile
//      );
//  }
//};



/// trampoline class for emf Filter
template<int D>
class PyFilter : public Filter<D>
{
  using Filter<D>::Filter;

  void solve( emf::Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Filter<D>,
      solve,
      tile
      );
  }
};



void bind_emf(py::module& m_sub)
{
    
  py::class_<
    emf::Grids,
    std::shared_ptr<emf::Grids>
            >(m_sub, "Grids")
    .def(py::init<int, int, int>())
    .def_readwrite("ex",   &emf::Grids::ex , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("ey",   &emf::Grids::ey , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("ez",   &emf::Grids::ez , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("bx",   &emf::Grids::bx , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("by",   &emf::Grids::by , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("bz",   &emf::Grids::bz , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jx",   &emf::Grids::jx , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jy",   &emf::Grids::jy , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jz",   &emf::Grids::jz , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("rho",  &emf::Grids::rho, py::return_value_policy::reference, py::keep_alive<1,0>());



  //--------------------------------------------------
  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  py::module m_3d = m_sub.def_submodule("threeD","3D specializations");

  /// General class for handling Maxwell's equations
  auto t1 = declare_tile<1>(m_1d, "Tile");
  auto t2 = declare_tile<2>(m_2d, "Tile");
  auto t3 = declare_tile<3>(m_3d, "Tile"); // defined below


  // FIXME extra debug additions/tests
  //t3.def_property("gs", 
  //  &emf::Tile<3>::get_grids2,
  //  &emf::Tile<3>::set_grids,
  //  py::return_value_policy::reference_internal, 
  //  py::keep_alive<0,1>());
  //t3.def("get_grids_ptr", &emf::Tile<3>::get_grids_ptr
  //    );
  //     // py::return_value_policy::reference_internal);


  // FIXME
  // Declare manually instead because there are too many differences
  //py::class_<emf::Tile<3>, corgi::Tile<3>, 
  //           std::shared_ptr<emf::Tile<3>>
  //          >(m_3d, "Tile")
  //  .def(py::init<size_t, size_t, size_t>())
  //  .def_readwrite("dx",         &emf::Tile<3>::dx)
  //  .def_readwrite("cfl",        &emf::Tile<3>::cfl)
  //  //.def_readwrite("gs",        &emf::Tile<3>::gs,
  //  //    py::return_value_policy::reference_internal, 
  //  //    py::keep_alive<0,1>())
  //  .def_property("gs", 
  //      &emf::Tile<3>::get_grids2,
  //      &emf::Tile<3>::set_grids,
  //      py::return_value_policy::reference_internal, 
  //      py::keep_alive<0,1>())
  //  .def("cycle_grids",            &emf::Tile<3>::cycle_grids)
  //  .def("clear_current",        &emf::Tile<3>::clear_current)
  //  .def("deposit_current",      &emf::Tile<3>::deposit_current)
  //  .def("update_boundaries",    &emf::Tile<3>::update_boundaries)
  //  .def("exchange_currents",    &emf::Tile<3>::exchange_currents)
  //  .def("get_grids_ptr",           &emf::Tile<3>::get_grids_ptr,
  //      py::return_value_policy::reference_internal)
  //  .def("get_grids",              &emf::Tile<3>::get_grids, 
  //      py::arg("i")=0,
  //      py::return_value_policy::reference,
  //      // keep alive for the lifetime of the grid
  //      //
  //      // pybind11:
  //      // argument indices start at one, while zero refers to the return 
  //      // value. For methods, index one refers to the implicit this pointer, 
  //      // while regular arguments begin at index two. 
  //      // py::keep_alive<nurse,patient>()
  //      py::keep_alive<1,0>()
  //      )


  // FIXME
  // TODO: testing grid-based tile generator instead
  //
  // example:
  //m.def("make_bell", []() { return new Bell; }, return_value_policy::reference);
  //m_3d.def("make_and_add_tile", [](
  //      corgi::Grid<3>& grid,
  //      int nx, int ny, int nz,
  //      corgi::internals::tuple_of<3, size_t> indices
  //      ) {

  //      //auto p = new emf::Tile<3>(nx,ny,nz);
  //      //std::shared_ptr<emf::Tile<3>> sp(p);
  //      //sp->index = indices;
  //      //grid.add_tile(sp, indices);
  //        
  //      std::shared_ptr<emf::Tile<3>> sp(new emf::Tile<3>(nx,ny,nz));
  //      grid.add_tile(sp, indices);

  //      //emf::Tile<3> ti(nx,ny,nz);
  //      //std::shared_ptr<emf::Tile<3>> sp(&ti);
  //      //grid.add_tile(sp, indices);

  //      return sp;

  //      //auto t = grid.get_tileptr(indices);
  //      //return t;
  //    },
  //      py::return_value_policy::reference,
  //      // keep alive for the lifetime of the grid
  //      //
  //      // pybind11:
  //      // argument indices start at one, while zero refers to the return 
  //      // value. For methods, index one refers to the implicit this pointer, 
  //      // while regular arguments begin at index two. 
  //      // py::keep_alive<nurse,patient>()
  //      py::keep_alive<1,0>()
  //    ); 


  //--------------------------------------------------
  // damping tiles

  auto td1_m1 = declare_TileDamped<1, -1>(m_1d, "TileDamped1D_LX");
  auto td1_p1 = declare_TileDamped<1, +1>(m_1d, "TileDamped1D_RX");

  auto td2_m1 = declare_TileDamped<2, -1>(m_2d, "TileDamped2D_LX");
  auto td2_p1 = declare_TileDamped<2, +1>(m_2d, "TileDamped2D_RX");
  auto td2_m2 = declare_TileDamped<2, -2>(m_2d, "TileDamped2D_LY");
  auto td2_p2 = declare_TileDamped<2, +2>(m_2d, "TileDamped2D_RY");

  auto td3_m1 = declare_TileDamped<3, -1>(m_3d, "TileDamped3D_LX");
  auto td3_p1 = declare_TileDamped<3, +1>(m_3d, "TileDamped3D_RX");
  auto td3_m2 = declare_TileDamped<3, -2>(m_3d, "TileDamped3D_LY");
  auto td3_p2 = declare_TileDamped<3, +2>(m_3d, "TileDamped3D_RY");
  auto td3_m3 = declare_TileDamped<3, -3>(m_3d, "TileDamped3D_LZ");
  auto td3_p3 = declare_TileDamped<3, +3>(m_3d, "TileDamped3D_RZ");

  //--------------------------------------------------
  // 1D Propagator bindings
  py::class_< emf::Propagator<1>, PyPropagator<1> > emfpropag1d(m_1d, "Propagator");
  emfpropag1d
    .def(py::init<>())
    .def("push_e",      &emf::Propagator<1>::push_e)
    .def("push_half_b", &emf::Propagator<1>::push_half_b);

  // fdtd2 propagator
  py::class_<emf::FDTD2<1>, Propagator<1>, PyFDTD2<1>>(m_1d, "FDTD2")
    .def(py::init<>())
    .def_readwrite("corr",     &emf::FDTD2<1>::corr);

  // fdtd4 propagator
  py::class_<emf::FDTD4<1>, Propagator<1>, PyFDTD4<1> >(m_1d, "FDTD4")
    .def_readwrite("corr",     &emf::FDTD4<1>::corr)
    .def(py::init<>());

  // fdtd2 propagator with perfectly matched ouer layer
  py::class_<emf::FDTD2_pml<1>> pml1d(m_1d, "FDTD2_pml", emfpropag1d);
  pml1d
    .def(py::init<>())
    .def_readwrite("cenx",     &emf::FDTD2_pml<1>::cenx)
    .def_readwrite("ceny",     &emf::FDTD2_pml<1>::ceny)
    .def_readwrite("cenz",     &emf::FDTD2_pml<1>::cenz)
    .def_readwrite("radx",     &emf::FDTD2_pml<1>::radx)
    .def_readwrite("rady",     &emf::FDTD2_pml<1>::rady)
    .def_readwrite("radz",     &emf::FDTD2_pml<1>::radz)
    .def_readwrite("rad_lim",  &emf::FDTD2_pml<1>::rad_lim)
    .def_readwrite("norm_abs", &emf::FDTD2_pml<1>::norm_abs)
    .def_readwrite("corr",     &emf::FDTD2_pml<1>::corr)
    .def_readwrite("mode",     &emf::FDTD2_pml<1>::mode)
    .def("push_e",             &emf::FDTD2_pml<1>::push_e)
    .def("push_half_b",        &emf::FDTD2_pml<1>::push_half_b);


  //--------------------------------------------------
  // 2D Propagator bindings
  py::class_< emf::Propagator<2>, PyPropagator<2> > emfpropag2d(m_2d, "Propagator");
  emfpropag2d
    .def(py::init<>())
    .def_readwrite("dt",&emf::Propagator<2>::dt)
    .def("push_e",      &emf::Propagator<2>::push_e)
    .def("push_half_b", &emf::Propagator<2>::push_half_b);

  // fdtd2 propagator
  py::class_<emf::FDTD2<2>>(m_2d, "FDTD2", emfpropag2d)
    .def_readwrite("corr",     &emf::FDTD2<2>::corr)
    .def(py::init<>());
    
  // fdtd2 propagator with perfectly matched ouer layer
  py::class_<emf::FDTD2_pml<2>> pml2d(m_2d, "FDTD2_pml", emfpropag2d);
  pml2d
    .def(py::init<>())
    .def_readwrite("cenx",     &emf::FDTD2_pml<2>::cenx)
    .def_readwrite("ceny",     &emf::FDTD2_pml<2>::ceny)
    .def_readwrite("cenz",     &emf::FDTD2_pml<2>::cenz)
    .def_readwrite("radx",     &emf::FDTD2_pml<2>::radx)
    .def_readwrite("rady",     &emf::FDTD2_pml<2>::rady)
    .def_readwrite("radz",     &emf::FDTD2_pml<2>::radz)
    .def_readwrite("rad_lim",  &emf::FDTD2_pml<2>::rad_lim)
    .def_readwrite("norm_abs", &emf::FDTD2_pml<2>::norm_abs)
    .def_readwrite("corr",     &emf::FDTD2_pml<2>::corr)
    .def_readwrite("mode",     &emf::FDTD2_pml<2>::mode)
    .def("push_e",             &emf::FDTD2_pml<2>::push_e)
    .def("push_half_b",        &emf::FDTD2_pml<2>::push_half_b);
    //.def("push_eb",            &emf::FDTD2_pml<2>::push_eb); // TODO not implemented


  // fdtd4 propagator
  py::class_<emf::FDTD4<2>, Propagator<2>, PyFDTD4<2> >(m_2d, "FDTD4")
    .def_readwrite("corr",     &emf::FDTD4<2>::corr)
    .def(py::init<>());


  //--------------------------------------------------
  // 3D Propagator bindings
  py::class_< emf::Propagator<3>, PyPropagator<3> > emfpropag3d(m_3d, "Propagator");
  emfpropag3d
    .def(py::init<>())
    .def("push_e",      &emf::Propagator<3>::push_e)
    .def("push_half_b", &emf::Propagator<3>::push_half_b);

  // fdtd2 propagator
  py::class_<emf::FDTD2<3>>(m_3d, "FDTD2", emfpropag3d)
    .def(py::init<>())
    .def_readwrite("corr",     &emf::FDTD2<3>::corr);

  // fdtd2 propagator with perfectly matched ouer layer
  py::class_<emf::FDTD2_pml<3>> pml3d(m_3d, "FDTD2_pml", emfpropag3d);
  pml3d
    .def(py::init<>())
    .def_readwrite("cenx",     &emf::FDTD2_pml<3>::cenx)
    .def_readwrite("ceny",     &emf::FDTD2_pml<3>::ceny)
    .def_readwrite("cenz",     &emf::FDTD2_pml<3>::cenz)
    .def_readwrite("radx",     &emf::FDTD2_pml<3>::radx)
    .def_readwrite("rady",     &emf::FDTD2_pml<3>::rady)
    .def_readwrite("radz",     &emf::FDTD2_pml<3>::radz)
    .def_readwrite("rad_lim",  &emf::FDTD2_pml<3>::rad_lim)
    .def_readwrite("norm_abs", &emf::FDTD2_pml<3>::norm_abs)
    .def_readwrite("corr",     &emf::FDTD2_pml<3>::corr)
    .def_readwrite("mode",     &emf::FDTD2_pml<3>::mode)
    .def("push_e",             &emf::FDTD2_pml<3>::push_e)
    .def("push_half_b",        &emf::FDTD2_pml<3>::push_half_b)
    .def("push_eb",            &emf::FDTD2_pml<3>::push_eb);


  // fdtd4 propagator
  py::class_<emf::FDTD4<3>, Propagator<3>, PyFDTD4<3> >(m_3d, "FDTD4")
    .def_readwrite("corr",     &emf::FDTD4<3>::corr)
    .def(py::init<>());


  // fdtd general propagator
  //py::class_<emf::FDTDGen<3>, Propagator<3>, PyFDTDGen<3> >(m_3d, "FDTDGen")
  py::class_<emf::FDTDGen<3>>(m_3d, "FDTDGen")
    .def_readwrite("corr",     &emf::FDTDGen<3>::corr)
    .def_readwrite("CXs",      &emf::FDTDGen<3>::CXs, py::return_value_policy::reference,py::keep_alive<1,0>())
    .def_readwrite("CYs",      &emf::FDTDGen<3>::CYs, py::return_value_policy::reference,py::keep_alive<1,0>())
    .def_readwrite("CZs",      &emf::FDTDGen<3>::CZs, py::return_value_policy::reference,py::keep_alive<1,0>())
    .def("push_e",             &emf::FDTDGen<3>::push_e)
    .def("push_half_b",        &emf::FDTDGen<3>::push_half_b)
    .def(py::init<>());


  //--------------------------------------------------
  // 1D, 2D, and 3D FILTERS

  // 1D filters
  py::class_< emf::Filter<1>, PyFilter<1> > emffilter1d(m_1d, "Filter");
  emffilter1d
    .def(py::init<int, int, int>())
    .def("solve", &emf::Filter<1>::solve);

  // digital filter
  py::class_<emf::Binomial2<1>>(m_1d, "Binomial2", emffilter1d)
    .def(py::init<int, int, int>())
    .def("solve",      &emf::Binomial2<1>::solve);
  
  // 2D Filter bindings
  py::class_< emf::Filter<2>, PyFilter<2> > emffilter2d(m_2d, "Filter");
  emffilter2d
    .def(py::init<int, int, int>())
    .def("solve", &emf::Filter<2>::solve);

  // digital filter
  // TODO: remove hack where we explicitly define solve (instead of use trampoline class)
  // overwriting the solve function from trampoline does not work atm for some weird reason.
  py::class_<emf::Binomial2<2>>(m_2d, "Binomial2", emffilter2d)
    .def(py::init<int, int, int>())
    .def("solve",      &emf::Binomial2<2>::solve);

  py::class_<emf::General3p<2>>(m_2d, "General3p", emffilter2d)
    .def(py::init<int, int, int>())
    .def_readwrite("alpha",    &emf::General3p<2>::alpha)
    .def("solve",              &emf::General3p<2>::solve);

  py::class_<emf::General3pStrided<2>>(m_2d, "General3pStrided", emffilter2d)
    .def(py::init<int, int, int>())
    .def_readwrite("alpha",    &emf::General3pStrided<2>::alpha)
    .def_readwrite("stride",   &emf::General3pStrided<2>::stride)
    .def("solve",              &emf::General3pStrided<2>::solve);


  py::class_<emf::Binomial2Strided2<2>>(m_2d, "Binomial2Strided2", emffilter2d)
    .def(py::init<int, int, int>())
    .def("solve",              &emf::Binomial2Strided2<2>::solve);

  py::class_<emf::Compensator2<2>>(m_2d, "Compensator2", emffilter2d)
    .def(py::init<int, int, int>())
    .def("solve",              &emf::Compensator2<2>::solve);



  // 3D filters
  py::class_< emf::Filter<3>, PyFilter<3> > emffilter3d(m_3d, "Filter");
  emffilter3d
    .def(py::init<int, int, int>())
    .def("solve", &emf::Filter<3>::solve);

  // digital filter
  py::class_<emf::Binomial2<3>>(m_3d, "Binomial2", emffilter3d)
    .def(py::init<int, int, int>())
    .def("solve",      &emf::Binomial2<3>::solve);



  //--------------------------------------------------
  // EM boundary conditions
    
  // 1D rotating conductor
  py::class_<emf::Conductor<1>>(m_1d, "Conductor")
    .def(py::init<>())
    .def_readwrite("B0",       &emf::Conductor<1>::B0)
    .def_readwrite("radius",   &emf::Conductor<1>::radius)
    .def_readwrite("period",   &emf::Conductor<1>::period)
    .def_readwrite("chi_om",   &emf::Conductor<1>::chi_om)
    .def_readwrite("chi_mu",   &emf::Conductor<1>::chi_mu)
    .def_readwrite("phase_mu", &emf::Conductor<1>::phase_mu)
    .def_readwrite("phase_om", &emf::Conductor<1>::phase_om)
    .def_readwrite("cenx",     &emf::Conductor<1>::cenx)
    .def_readwrite("ceny",     &emf::Conductor<1>::ceny)
    .def_readwrite("cenz",     &emf::Conductor<1>::cenz)
    .def_readwrite("delta",    &emf::Conductor<1>::delta)
    .def_readwrite("radius_pc",&emf::Conductor<1>::radius_pc)
    .def_readwrite("delta_pc", &emf::Conductor<1>::delta_pc)
    .def_readwrite("Nx",       &emf::Conductor<1>::Nx)
    .def_readwrite("Ny",       &emf::Conductor<1>::Ny)
    .def_readwrite("Nz",       &emf::Conductor<1>::Nz)
    .def("insert_em",          &emf::Conductor<1>::insert_em)
    .def("update_b",           &emf::Conductor<1>::update_b)
    .def("update_j",           &emf::Conductor<1>::update_j)
    .def("update_e",           &emf::Conductor<1>::update_e);

    
  // 2D rotating conductor
  py::class_<emf::Conductor<2>>(m_2d, "Conductor")
    .def(py::init<>())
    .def_readwrite("B0",       &emf::Conductor<2>::B0)
    .def_readwrite("radius",   &emf::Conductor<2>::radius)
    .def_readwrite("period",   &emf::Conductor<2>::period)
    .def_readwrite("chi_om",   &emf::Conductor<2>::chi_om)
    .def_readwrite("chi_mu",   &emf::Conductor<2>::chi_mu)
    .def_readwrite("phase_mu", &emf::Conductor<2>::phase_mu)
    .def_readwrite("phase_om", &emf::Conductor<2>::phase_om)
    .def_readwrite("cenx",     &emf::Conductor<2>::cenx)
    .def_readwrite("ceny",     &emf::Conductor<2>::ceny)
    .def_readwrite("cenz",     &emf::Conductor<2>::cenz)
    .def_readwrite("delta",    &emf::Conductor<2>::delta)
    .def_readwrite("radius_pc",&emf::Conductor<2>::radius_pc)
    .def_readwrite("delta_pc", &emf::Conductor<2>::delta_pc)
    .def_readwrite("Nx",       &emf::Conductor<2>::Nx)
    .def_readwrite("Ny",       &emf::Conductor<2>::Ny)
    .def_readwrite("Nz",       &emf::Conductor<2>::Nz)
    .def("insert_em",          &emf::Conductor<2>::insert_em)
    .def("update_b",           &emf::Conductor<2>::update_b)
    .def("update_j",           &emf::Conductor<2>::update_j)
    .def("update_e",           &emf::Conductor<2>::update_e);


  // 3D rotating conductor
  py::class_<emf::Conductor<3>>(m_3d, "Conductor")
    .def(py::init<>())
    .def_readwrite("B0",       &emf::Conductor<3>::B0)
    .def_readwrite("radius",   &emf::Conductor<3>::radius)
    .def_readwrite("period",   &emf::Conductor<3>::period)
    .def_readwrite("chi_om",   &emf::Conductor<3>::chi_om)
    .def_readwrite("chi_mu",   &emf::Conductor<3>::chi_mu)
    .def_readwrite("phase_mu", &emf::Conductor<3>::phase_mu)
    .def_readwrite("phase_om", &emf::Conductor<3>::phase_om)
    .def_readwrite("cenx",     &emf::Conductor<3>::cenx)
    .def_readwrite("ceny",     &emf::Conductor<3>::ceny)
    .def_readwrite("cenz",     &emf::Conductor<3>::cenz)
    .def_readwrite("delta",    &emf::Conductor<3>::delta)
    .def_readwrite("radius_pc",&emf::Conductor<3>::radius_pc)
    .def_readwrite("delta_pc", &emf::Conductor<3>::delta_pc)
    .def_readwrite("Nx",       &emf::Conductor<3>::Nx)
    .def_readwrite("Ny",       &emf::Conductor<3>::Ny)
    .def_readwrite("Nz",       &emf::Conductor<3>::Nz)
    .def("insert_em",          &emf::Conductor<3>::insert_em)
    .def("update_e",           &emf::Conductor<3>::update_e)
    .def("update_j",           &emf::Conductor<3>::update_j)
    .def("update_b",           &emf::Conductor<3>::update_b);


  //--------------------------------------------------
  // Snapshot IO 
  
  // 1D 
  py::class_<h5io::FieldsWriter<1>>(m_1d, "FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldsWriter<1>::write) 
    .def("get_slice", [](h5io::FieldsWriter<1> &s, int k)
            {
                const auto nx = static_cast<pybind11::ssize_t>( s.nx );
                const auto ny = static_cast<pybind11::ssize_t>( s.ny );
                auto v = pybind11::array_t<float>( {nx, ny}, s.arrs[k].data() );
                return v;
            });

  // 2D 
  py::class_<h5io::FieldsWriter<2>>(m_2d, "FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldsWriter<2>::write) 
    .def("get_slice", [](h5io::FieldsWriter<2> &s, int k)
            {
                //const auto N = static_cast<pybind11::ssize_t>(s.arrs[k].size());
                const auto nx = static_cast<pybind11::ssize_t>( s.nx );
                const auto ny = static_cast<pybind11::ssize_t>( s.ny );
                auto v = pybind11::array_t<float>( {nx, ny}, s.arrs[k].data() );
                return v;
            });

  // 3D 
  py::class_<h5io::FieldsWriter<3>>(m_3d, "FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldsWriter<3>::write);

  // 3D; root only field storage
  py::class_<h5io::MasterFieldsWriter<3>>(m_3d, "MasterFieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::MasterFieldsWriter<3>::write);

  // slice writer; only in 3D
  py::class_<h5io::FieldSliceWriter>(m_3d, "FieldSliceWriter")
    .def_readwrite("ind",  &h5io::FieldSliceWriter::ind)
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",        &h5io::FieldSliceWriter::write)
    .def("get_slice", [](h5io::FieldSliceWriter &s, int k)
            {
                //const auto N = static_cast<pybind11::ssize_t>(s.arrs[k].size());
                const auto nx = static_cast<pybind11::ssize_t>( s.nx );
                const auto ny = static_cast<pybind11::ssize_t>( s.ny );
                auto v = pybind11::array_t<float>( {nx, ny}, s.arrs[k].data() );
                return v;
            });

  //--------------------------------------------------
  // Full IO 

  // 1D
  m_1d.def("read_grids",        &emf::read_grids<1>);
  m_1d.def("write_grids",       &emf::write_grids<1>);


  // 2D
  m_2d.def("write_grids",        &emf::write_grids<2>);
  m_2d.def("read_grids",         &emf::read_grids<2>);


  // 3D
  m_3d.def("write_grids",        &emf::write_grids<3>);
  m_3d.def("read_grids",         &emf::read_grids<3>);



}

} // end of namespace emf
