
#include "py_submodules.h"


#include "../definitions.h"
#include "../tools/mesh.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

//--------------------------------------------------
  
namespace fields{
  
// generator for damped tile for various directions
template<size_t D>
auto declare_Tile(
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
    .def("cycleYee",           &fields::Tile<D>::cycleYee)
    .def("cycleCurrent",       &fields::Tile<D>::cycleCurrent)
    .def("pushE",              &fields::Tile<D>::pushE)
    .def("pushHalfB",          &fields::Tile<D>::pushHalfB)
    .def("depositCurrent",     &fields::Tile<D>::depositCurrent)
    .def("getYee",             &fields::Tile<D>::getYee, py::return_value_policy::reference)
    .def("getAnalysis",        &fields::Tile<D>::getAnalysis, py::return_value_policy::reference)
    .def("addAnalysisSpecies", &fields::Tile<D>::addAnalysisSpecies)
    .def("updateBoundaries",   &fields::Tile<D>::updateBoundaries)
    .def("exchangeCurrents",   &fields::Tile<D>::exchangeCurrents);
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
             corgi::Tile<D>,
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
  .def("dampFields",         &fields::damping::Tile<D,S>::dampFields);
}



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
    .def_readwrite("rho",  &fields::YeeLattice::rho);

  //--------------------------------------------------

  py::class_<fields::PlasmaMomentLattice>(m_sub, "PlasmaMomentLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("rho",      &fields::PlasmaMomentLattice::rho)
    .def_readwrite("mgamma",   &fields::PlasmaMomentLattice::mgamma)
    .def_readwrite("Vx",       &fields::PlasmaMomentLattice::Vx)
    .def_readwrite("Vy",       &fields::PlasmaMomentLattice::Vy)
    .def_readwrite("Vz",       &fields::PlasmaMomentLattice::Vz)
    .def_readwrite("Tx",       &fields::PlasmaMomentLattice::Tx)
    .def_readwrite("Ty",       &fields::PlasmaMomentLattice::Ty)
    .def_readwrite("Tz",       &fields::PlasmaMomentLattice::Tz)
    .def_readwrite("ekin",     &fields::PlasmaMomentLattice::ekin);


  //--------------------------------------------------
  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");

  /// General class for handling Maxwell's equations
  auto t1 = declare_Tile<1>(m_1d, "Tile");
  auto t2 = declare_Tile<2>(m_2d, "Tile");
  //auto t3 = declare_Tile<3>(m, "Tile");


  //--------------------------------------------------
  // damping tiles

  auto td1_m1 = declare_TileDamped<1, -1>(m_1d, "TileDamped1D_LX");
  auto td1_p1 = declare_TileDamped<1, +1>(m_1d, "TileDamped1D_RX");

  auto td2_m1 = declare_TileDamped<2, -1>(m_2d, "TileDamped2D_LX");
  auto td2_p1 = declare_TileDamped<2, +1>(m_2d, "TileDamped2D_RX");
  auto td2_m2 = declare_TileDamped<2, -2>(m_2d, "TileDamped2D_LY");
  auto td2_p2 = declare_TileDamped<2, +2>(m_2d, "TileDamped2D_RY");


}

} // end of namespace fields
