#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#include "../ffe/tile.h"

#include "../ffe/currents/current.h"
#include "../ffe/currents/driftcurrent.h"
#include "../ffe/skinny_yee.h"


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
    .def("get_step", [](ffe::Tile<D>& s, int n)
        -> SkinnyYeeLattice& 
        {
        if(n < 4) return s.get_step(n);

        // else
        throw py::index_error();
        },
          py::return_value_policy::reference
    );

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

  // skinny version of the Yee lattice with only (e and b meshes)
  py::class_<ffe::SkinnyYeeLattice>(m_sub, "SkinnyYeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex",   &ffe::SkinnyYeeLattice::ex)
    .def_readwrite("ey",   &ffe::SkinnyYeeLattice::ey)
    .def_readwrite("ez",   &ffe::SkinnyYeeLattice::ez)
    .def_readwrite("bx",   &ffe::SkinnyYeeLattice::bx)
    .def_readwrite("by",   &ffe::SkinnyYeeLattice::by)
    .def_readwrite("bz",   &ffe::SkinnyYeeLattice::bz)
    .def("set_yee",        &ffe::SkinnyYeeLattice::set_yee)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self +  py::self)
    .def(py::self -  py::self)
    .def(py::self *  float())
    .def(py::self /  float());


  m_sub.def("set_step", [](fields::YeeLattice& yee, ffe::SkinnyYeeLattice skyee)
      -> void 
      {
        yee.ex = skyee.ex;
        yee.ey = skyee.ey;
        yee.ez = skyee.ez;
        
        yee.bx = skyee.bx;
        yee.by = skyee.by;
        yee.bz = skyee.bz;
        return;
      }
  );


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
