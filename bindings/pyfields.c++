#include <string>

#include "py_submodules.h"

#include "../corgi/internals.h"

#include "../definitions.h"
#include "../tools/mesh.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

#include "../em-fields/propagator/propagator.h"
#include "../em-fields/propagator/fdtd2.h"
#include "../em-fields/propagator/fdtd4.h"

#include "../em-fields/filters/filter.h"
#include "../em-fields/filters/digital.h"

#include "../io/writers/writer.h"
#include "../io/writers/fields.h"
#include "../io/snapshots/fields.h"
#include "../io/snapshots/field_slices.h"
#include "../io/tasker.h"


//--------------------------------------------------
  
namespace fields{

  namespace py = pybind11;

  
// generator for damped tile for various directions
template<int D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{
  return py::class_<
             fields::Tile<D>,
              corgi::Tile<D>, 
              std::shared_ptr<fields::Tile<D>>
            >(m, pyclass_name.c_str() )
    .def(py::init<int, int, int>())
    .def_readwrite("cfl",      &fields::Tile<D>::cfl)
    .def("cycle_yee",           &fields::Tile<D>::cycle_yee)
    .def("clear_current",       &fields::Tile<D>::clear_current)
    .def("deposit_current",     &fields::Tile<D>::deposit_current)
    .def("update_boundaries",   &fields::Tile<D>::update_boundaries)
    .def("exchange_currents",   &fields::Tile<D>::exchange_currents)

    .def("get_yee",             &fields::Tile<D>::get_yee,
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
  // using Class = fields::TileDamped<D,S>; 
  // does not function properly; maybe not triggering template?
  // have to use explicit name instead like this

  return py::class_<
             fields::damping::Tile<D,S>,
             fields::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<fields::damping::Tile<D,S>>
          >(m, 
            pyclass_name.c_str(),
            py::multiple_inheritance()
            )
  .def(py::init<int, int, int>())
  .def_readwrite("ex_ref",   &fields::damping::Tile<D,S>::ex_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("ey_ref",   &fields::damping::Tile<D,S>::ey_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("ez_ref",   &fields::damping::Tile<D,S>::ez_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("bx_ref",   &fields::damping::Tile<D,S>::bx_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("by_ref",   &fields::damping::Tile<D,S>::by_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("bz_ref",   &fields::damping::Tile<D,S>::bz_ref, py::return_value_policy::reference,py::keep_alive<1,0>())
  .def_readwrite("fld1",     &fields::damping::Tile<D,S>::fld1)
  .def_readwrite("fld2",     &fields::damping::Tile<D,S>::fld2)
  .def_readwrite("ksupp",    &fields::damping::Tile<D,S>::ksupp)
  .def("damp_fields",        &fields::damping::Tile<D,S>::damp_fields);
}


/// trampoline class for fields Propagator
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
template<int D>
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
    
  py::class_<
    fields::YeeLattice,
    std::shared_ptr<fields::YeeLattice>
            >(m_sub, "YeeLattice")
    .def(py::init<int, int, int>())
    .def_readwrite("ex",   &fields::YeeLattice::ex , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("ey",   &fields::YeeLattice::ey , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("ez",   &fields::YeeLattice::ez , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("bx",   &fields::YeeLattice::bx , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("by",   &fields::YeeLattice::by , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("bz",   &fields::YeeLattice::bz , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jx",   &fields::YeeLattice::jx , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jy",   &fields::YeeLattice::jy , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("jz",   &fields::YeeLattice::jz , py::return_value_policy::reference, py::keep_alive<1,0>())
    .def_readwrite("rho",  &fields::YeeLattice::rho, py::return_value_policy::reference, py::keep_alive<1,0>());



  //--------------------------------------------------
  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  py::module m_3d = m_sub.def_submodule("threeD","3D specializations");

  /// General class for handling Maxwell's equations
  auto t1 = declare_tile<1>(m_1d, "Tile");
  auto t2 = declare_tile<2>(m_2d, "Tile");
  auto t3 = declare_tile<3>(m_3d, "Tile"); // defined below


  // FIXME extra debug additions/tests
  //t3.def_property("yee", 
  //  &fields::Tile<3>::get_yee2,
  //  &fields::Tile<3>::set_yee,
  //  py::return_value_policy::reference_internal, 
  //  py::keep_alive<0,1>());
  //t3.def("get_yeeptr", &fields::Tile<3>::get_yeeptr
  //    );
  //     // py::return_value_policy::reference_internal);


  // FIXME
  // Declare manually instead because there are too many differences
  //py::class_<fields::Tile<3>, corgi::Tile<3>, 
  //           std::shared_ptr<fields::Tile<3>>
  //          >(m_3d, "Tile")
  //  .def(py::init<size_t, size_t, size_t>())
  //  .def_readwrite("dx",         &fields::Tile<3>::dx)
  //  .def_readwrite("cfl",        &fields::Tile<3>::cfl)
  //  //.def_readwrite("yee",        &fields::Tile<3>::yee,
  //  //    py::return_value_policy::reference_internal, 
  //  //    py::keep_alive<0,1>())
  //  .def_property("yee", 
  //      &fields::Tile<3>::get_yee2,
  //      &fields::Tile<3>::set_yee,
  //      py::return_value_policy::reference_internal, 
  //      py::keep_alive<0,1>())
  //  .def("cycle_yee",            &fields::Tile<3>::cycle_yee)
  //  .def("clear_current",        &fields::Tile<3>::clear_current)
  //  .def("deposit_current",      &fields::Tile<3>::deposit_current)
  //  .def("update_boundaries",    &fields::Tile<3>::update_boundaries)
  //  .def("exchange_currents",    &fields::Tile<3>::exchange_currents)
  //  .def("get_yeeptr",           &fields::Tile<3>::get_yeeptr,
  //      py::return_value_policy::reference_internal)
  //  .def("get_yee",              &fields::Tile<3>::get_yee, 
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

  //      //auto p = new fields::Tile<3>(nx,ny,nz);
  //      //std::shared_ptr<fields::Tile<3>> sp(p);
  //      //sp->index = indices;
  //      //grid.add_tile(sp, indices);
  //        
  //      std::shared_ptr<fields::Tile<3>> sp(new fields::Tile<3>(nx,ny,nz));
  //      grid.add_tile(sp, indices);

  //      //fields::Tile<3> ti(nx,ny,nz);
  //      //std::shared_ptr<fields::Tile<3>> sp(&ti);
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
  // 3D Propagator bindings
  py::class_< fields::Propagator<3>, PyPropagator<3> > fieldspropag3d(m_3d, "Propagator");
  fieldspropag3d
    .def(py::init<>())
    .def("push_e",      &fields::Propagator<3>::push_e)
    .def("push_half_b", &fields::Propagator<3>::push_half_b);

  // fdtd2 propagator
  py::class_<fields::FDTD2<3>>(m_3d, "FDTD2", fieldspropag3d)
    .def(py::init<>())
    .def_readwrite("corr",     &fields::FDTD2<3>::corr);


  // fdtd4 propagator
  py::class_<fields::FDTD4<3>, Propagator<3>, PyFDTD4<3> >(m_3d, "FDTD4")
    .def_readwrite("corr",     &fields::FDTD4<3>::corr)
    .def(py::init<>());


  //--------------------------------------------------
  // 2D Filter bindings
  py::class_< fields::Filter<2>, PyFilter<2> > fieldsfilter2d(m_2d, "Filter");
  fieldsfilter2d
    .def(py::init<int, int, int>())
    .def("solve", &fields::Filter<2>::solve);

  // digital filter
  // TODO: remove hack where we explicitly define solve (instead of use trampoline class)
  // overwriting the solve function from trampoline does not work atm for some weird reason.
  py::class_<fields::Binomial2<2>>(m_2d, "Binomial2", fieldsfilter2d)
    .def(py::init<int, int, int>())
    .def("solve",      &fields::Binomial2<2>::solve);

  py::class_<fields::General3p<2>>(m_2d, "General3p", fieldsfilter2d)
    .def(py::init<int, int, int>())
    .def_readwrite("alpha",    &fields::General3p<2>::alpha)
    .def("solve",              &fields::General3p<2>::solve);

  py::class_<fields::General3pStrided<2>>(m_2d, "General3pStrided", fieldsfilter2d)
    .def(py::init<int, int, int>())
    .def_readwrite("alpha",    &fields::General3pStrided<2>::alpha)
    .def_readwrite("stride",   &fields::General3pStrided<2>::stride)
    .def("solve",              &fields::General3pStrided<2>::solve);


  py::class_<fields::Binomial2Strided2<2>>(m_2d, "Binomial2Strided2", fieldsfilter2d)
    .def(py::init<int, int, int>())
    .def("solve",              &fields::Binomial2Strided2<2>::solve);

  py::class_<fields::Compensator2<2>>(m_2d, "Compensator2", fieldsfilter2d)
    .def(py::init<int, int, int>())
    .def("solve",              &fields::Compensator2<2>::solve);


  // 3D filters

  py::class_< fields::Filter<3>, PyFilter<3> > fieldsfilter3d(m_3d, "Filter");
  fieldsfilter3d
    .def(py::init<int, int, int>())
    .def("solve", &fields::Filter<3>::solve);

  // digital filter
  py::class_<fields::Binomial2<3>>(m_3d, "Binomial2", fieldsfilter3d)
    .def(py::init<int, int, int>())
    .def("solve",      &fields::Binomial2<3>::solve);


  //--------------------------------------------------
  // Snapshot IO 

  // 2D 
  py::class_<h5io::FieldsWriter<2>>(m_2d, "FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldsWriter<2>::write);

  // 3D 
  py::class_<h5io::FieldsWriter<3>>(m_3d, "FieldsWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldsWriter<3>::write);


  // slice writer; only in 3D
  py::class_<h5io::FieldSliceWriter>(m_3d, "FieldSliceWriter")
    .def(py::init<const std::string&, int, int, int, int, int, int, int, int, int>())
    .def("write",   &h5io::FieldSliceWriter::write);


  //--------------------------------------------------
  // Full IO 

  // 1D
  m_1d.def("read_yee",        &fields::read_yee<1>);
  m_1d.def("write_yee",       &fields::write_yee<1>);


  // 2D
  m_2d.def("write_yee",        &fields::write_yee<2>);
  m_2d.def("read_yee",         &fields::read_yee<2>);


  // 3D
  m_3d.def("write_yee",        &fields::write_yee<3>);
  m_3d.def("read_yee",         &fields::read_yee<3>);



}

} // end of namespace fields
