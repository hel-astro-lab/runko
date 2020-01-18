#include <string>

#include "py_submodules.h"

#include "../definitions.h"
#include "../tools/mesh.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

#include "../em-fields/propagator/propagator.h"
#include "../em-fields/propagator/fdtd2.h"
#include "../em-fields/propagator/fdtd4.h"

#include "../em-fields/filters/filter.h"
#include "../em-fields/filters/digital.h"

#include "../io/quick_writer.h"

//--------------------------------------------------
  
namespace fields{

  namespace py = pybind11;

  
// generator for damped tile for various directions
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return py::class_<
             fields::Tile<D>,
              corgi::Tile<D>, 
             std::shared_ptr<fields::Tile<D> >
            >(m, pyclass_name.c_str() )
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("dx",       &fields::Tile<D>::dx)
    .def_readwrite("cfl",      &fields::Tile<D>::cfl)
    .def("cycle_yee",           &fields::Tile<D>::cycle_yee)
    .def("cycle_current",       &fields::Tile<D>::cycle_current)
    .def("clear_current",       &fields::Tile<D>::clear_current)
    .def("deposit_current",     &fields::Tile<D>::deposit_current)
    .def("add_analysis_species", &fields::Tile<D>::add_analysis_species)
    .def("update_boundaries",   &fields::Tile<D>::update_boundaries)
    .def("exchange_currents",   &fields::Tile<D>::exchange_currents)
    .def("get_yee",             &fields::Tile<D>::get_yee, 
                                py::arg("i")=0,
                                py::return_value_policy::reference)
    .def("get_analysis",        &fields::Tile<D>::get_analysis, 
                                py::arg("i")=0,
                                py::return_value_policy::reference);
        
}


// generator for damped tile for various directions
template<size_t D, int S>
auto declare_TileDamped(
    py::module& m,
    const std::string& pyclass_name) 
{
  //using Class = fields::TileDamped<D,S>; // does not function properly; maybe not triggering template?
  //have to use explicit name instead like this

  return py::class_<
             fields::damping::Tile<D,S>,
             fields::Tile<D>,
             std::shared_ptr<fields::damping::Tile<D,S>>
          >(m, pyclass_name.c_str() )
  .def(py::init<size_t, size_t, size_t>())
  .def_readwrite("ex_ref",   &fields::damping::Tile<D,S>::ex_ref, py::return_value_policy::reference)
  .def_readwrite("ey_ref",   &fields::damping::Tile<D,S>::ey_ref, py::return_value_policy::reference)
  .def_readwrite("ez_ref",   &fields::damping::Tile<D,S>::ez_ref, py::return_value_policy::reference)
  .def_readwrite("bx_ref",   &fields::damping::Tile<D,S>::bx_ref, py::return_value_policy::reference)
  .def_readwrite("by_ref",   &fields::damping::Tile<D,S>::by_ref, py::return_value_policy::reference)
  .def_readwrite("bz_ref",   &fields::damping::Tile<D,S>::bz_ref, py::return_value_policy::reference)
  .def_readwrite("fld1",     &fields::damping::Tile<D,S>::fld1)
  .def_readwrite("fld2",     &fields::damping::Tile<D,S>::fld2)
  .def("damp_fields",         &fields::damping::Tile<D,S>::damp_fields);
}


/// trampoline class for fields Propagator
template<size_t D>
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



/// trampoline class for fields Filter
template<size_t D>
class PyFilter : public Filter<D>
{
  using Filter<D>::Filter;

  void solve( fields::Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Filter<D>,
      solve,
      tile
      );
  }
};



void bind_fields(py::module& m_sub)
{
    
  //--------------------------------------------------

  py::class_<fields::YeeLattice>(m_sub, "YeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex",   &fields::YeeLattice::ex)
    .def_readwrite("ey",   &fields::YeeLattice::ey)
    .def_readwrite("ez",   &fields::YeeLattice::ez)
    .def_readwrite("bx",   &fields::YeeLattice::bx)
    .def_readwrite("by",   &fields::YeeLattice::by)
    .def_readwrite("bz",   &fields::YeeLattice::bz)
    .def_readwrite("jx",   &fields::YeeLattice::jx)
    .def_readwrite("jy",   &fields::YeeLattice::jy)
    .def_readwrite("jz",   &fields::YeeLattice::jz)
    .def_readwrite("jx1",  &fields::YeeLattice::jx1)
    .def_readwrite("jy1",  &fields::YeeLattice::jy1)
    .def_readwrite("jz1",  &fields::YeeLattice::jz1)
    .def_readwrite("rho",  &fields::YeeLattice::rho);

  //--------------------------------------------------

  py::class_<fields::PlasmaMomentLattice>(m_sub, "PlasmaMomentLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("rho",      &fields::PlasmaMomentLattice::rho)
    .def_readwrite("edens",    &fields::PlasmaMomentLattice::edens)
    .def_readwrite("temp",     &fields::PlasmaMomentLattice::temp)
    .def_readwrite("Vx",       &fields::PlasmaMomentLattice::Vx)
    .def_readwrite("Vy",       &fields::PlasmaMomentLattice::Vy)
    .def_readwrite("Vz",       &fields::PlasmaMomentLattice::Vz)
    .def_readwrite("momx",     &fields::PlasmaMomentLattice::momx)
    .def_readwrite("momy",     &fields::PlasmaMomentLattice::momy)
    .def_readwrite("momz",     &fields::PlasmaMomentLattice::momz)
    .def_readwrite("pressx",   &fields::PlasmaMomentLattice::pressx)
    .def_readwrite("pressy",   &fields::PlasmaMomentLattice::pressy)
    .def_readwrite("pressz",   &fields::PlasmaMomentLattice::pressz)
    .def_readwrite("shearxy",  &fields::PlasmaMomentLattice::shearxy)
    .def_readwrite("shearxz",  &fields::PlasmaMomentLattice::shearxz)
    .def_readwrite("shearyz",  &fields::PlasmaMomentLattice::shearyz);



  //--------------------------------------------------
  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");

  /// General class for handling Maxwell's equations
  auto t1 = declare_tile<1>(m_1d, "Tile");
  auto t2 = declare_tile<2>(m_2d, "Tile");
  //auto t3 = declare_tile<3>(m, "Tile");


  //--------------------------------------------------
  // damping tiles

  auto td1_m1 = declare_TileDamped<1, -1>(m_1d, "TileDamped1D_LX");
  auto td1_p1 = declare_TileDamped<1, +1>(m_1d, "TileDamped1D_RX");

  auto td2_m1 = declare_TileDamped<2, -1>(m_2d, "TileDamped2D_LX");
  auto td2_p1 = declare_TileDamped<2, +1>(m_2d, "TileDamped2D_RX");
  auto td2_m2 = declare_TileDamped<2, -2>(m_2d, "TileDamped2D_LY");
  auto td2_p2 = declare_TileDamped<2, +2>(m_2d, "TileDamped2D_RY");


  //--------------------------------------------------
  // 1D Propagator bindings
  py::class_< fields::Propagator<1>, PyPropagator<1> >(m_1d, "Propagator")
    .def(py::init<>())
    .def("push_e",      &fields::Propagator<1>::push_e)
    .def("push_half_b", &fields::Propagator<1>::push_half_b);

  // fdtd2 propagator
  py::class_<fields::FDTD2<1>, Propagator<1>, PyFDTD2<1>>(m_1d, "FDTD2")
    .def(py::init<>())
    .def_readwrite("corr",     &fields::FDTD2<1>::corr);


  //--------------------------------------------------
  // 2D Propagator bindings
  py::class_< fields::Propagator<2>, PyPropagator<2> >(m_2d, "Propagator")
    .def(py::init<>())
    .def_readwrite("dt",&fields::Propagator<2>::dt)
    .def("push_e",      &fields::Propagator<2>::push_e)
    .def("push_half_b", &fields::Propagator<2>::push_half_b);

  // fdtd2 propagator
  py::class_<fields::FDTD2<2>, Propagator<2>, PyFDTD2<2> >(m_2d, "FDTD2")
    .def_readwrite("corr",     &fields::FDTD2<2>::corr)
    .def(py::init<>());

  // fdtd4 propagator
  py::class_<fields::FDTD4<2>, Propagator<2>, PyFDTD4<2> >(m_2d, "FDTD4")
    .def_readwrite("corr",     &fields::FDTD4<2>::corr)
    .def(py::init<>());


  //--------------------------------------------------
  // 2D Filter bindings
  py::class_< fields::Filter<2>, PyFilter<2> > fieldsfilter2d(m_2d, "Filter");
  fieldsfilter2d
    .def(py::init<size_t, size_t, size_t>())
    .def("solve", &fields::Filter<2>::solve);

  // digital filter
  // TODO: remove hack where we explicitly define solve (instead of use trampoline class)
  // overwriting the solve function from trampoline does not work atm for some weird reason.
  py::class_<fields::Binomial2<2>>(m_2d, "Binomial2", fieldsfilter2d)
    .def(py::init<size_t, size_t, size_t>())
    .def("solve",      &fields::Binomial2<2>::solve);

  py::class_<fields::General3p<2>>(m_2d, "General3p", fieldsfilter2d)
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("alpha",    &fields::General3p<2>::alpha)
    .def("solve",              &fields::General3p<2>::solve);

  py::class_<fields::General3pStrided<2>>(m_2d, "General3pStrided", fieldsfilter2d)
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("alpha",    &fields::General3pStrided<2>::alpha)
    .def_readwrite("stride",   &fields::General3pStrided<2>::stride)
    .def("solve",              &fields::General3pStrided<2>::solve);


  py::class_<fields::Binomial2Strided2<2>>(m_2d, "Binomial2Strided2", fieldsfilter2d)
    .def(py::init<size_t, size_t, size_t>())
    .def("solve",              &fields::Binomial2Strided2<2>::solve);

  py::class_<fields::Compensator2<2>>(m_2d, "Compensator2", fieldsfilter2d)
    .def(py::init<size_t, size_t, size_t>())
    .def("solve",              &fields::Compensator2<2>::solve);



  //--------------------------------------------------
  // Quick IO 

  py::class_<h5io::QuickWriter<2>>(m_2d, "QuickWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::QuickWriter<2>::write);


}

} // end of namespace fields
