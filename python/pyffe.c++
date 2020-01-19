#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "../ffe/tile.h"

#include "../ffe/currents/current.h"
#include "../ffe/currents/driftcurrent.h"


namespace ffe {

//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<ffe::Tile<D>, 
             fields::Tile<D>,
             std::shared_ptr<ffe::Tile<D>>
             >(m, pyclass_name.c_str())
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("dx",        &ffe::Tile<D>::dx)
    .def_readwrite("cfl",       &ffe::Tile<D>::cfl)
    .def("compute_perp_current",    &ffe::Tile<D>::compute_perp_current)
    .def("subtract_parallel_e",     &ffe::Tile<D>::subtract_parallel_e);
}

//--------------------------------------------------
//
/// trampoline class for Current solver
template<size_t D>
class PyCurrent : public Current<D>
{
  using Current<D>::Current;

  void comp_drift_cur( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Current<D>,
      comp_drift_cur,
      tile
      );
  }

  void comp_parallel_cur( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Current<D>,
      comp_parallel_cur,
      tile
      );
  }

  void limiter( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      Current<D>,
      limiter,
      tile
      );
  }

};


template<size_t D>
class PyDriftCurrent : public DriftCurrent<D>
{
  using DriftCurrent<D>::DriftCurrent;

  void comp_drift_cur( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      DriftCurrent<D>,
      comp_drift_cur,
      tile
      );
  }

  void comp_parallel_cur( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      DriftCurrent<D>,
      comp_parallel_cur,
      tile
      );
  }

  void limiter( Tile<D>& tile ) override {
  PYBIND11_OVERLOAD_PURE(
      void,
      DriftCurrent<D>,
      limiter,
      tile
      );
  }
};




// python bindings for plasma classes & functions
void bind_ffe(py::module& m_sub)
{
  //--------------------------------------------------
  // 1D bindings
  //py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");

  //--------------------------------------------------
  // 2D bindings
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  auto t2 = ffe::declare_tile<2>(m_2d, "Tile");


  //--------------------------------------------------
  // 3D bindings
  //py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");
  //auto t3 = ffe::declare_tile<3>(m_2d, "Tile");


  //--------------------------------------------------
  // 2D Current solver bindings
  py::class_< ffe::Current<2>, PyCurrent<2> > currentcalc2d(m_2d, "Current");
  currentcalc2d
    .def(py::init<int, int, int>())
    .def("comp_drift_cur",      &ffe::Current<2>::comp_drift_cur)
    .def("comp_parallel_cur",   &ffe::Current<2>::comp_parallel_cur)
    .def("limiter",             &ffe::Current<2>::limiter);

  // Drift current solver
  py::class_<ffe::DriftCurrent<2>, ffe::Current<2>, PyDriftCurrent<2> >(m_2d, "DriftCurrent")
    .def(py::init<int,int,int>());



}

} // end of ns ffe
