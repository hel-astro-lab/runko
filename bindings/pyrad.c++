#include "py_submodules.h"

#include "../qed/photon.h"
#include "../qed/tile.h"


namespace qed {

//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<qed::Tile<D>,
             pic::Tile<D>, 
             fields::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<qed::Tile<D>>
             >(m, 
               pyclass_name.c_str(),
               py::multiple_inheritance()
               )
    .def(py::init<int, int, int>())
    .def("get_bucket",       &qed::Tile<D>::get_bucket, 
        py::return_value_policy::reference)
    .def("push_back",       &qed::Tile<D>::push_back);
}



//--------------------------------------------------

// python bindings for radiation classes & functions
void bind_qed(py::module& m_sub)
{

  py::class_<qed::PhotonContainer, pic::ParticleContainer<3>>(m_sub, "PhotonContainer")
    .def(py::init<>())
    .def("add_particle",  (void (qed::PhotonContainer::*)
          ( std::vector<float_p>, std::vector<float_p>, float_p, float_p ) ) 
            &qed::PhotonContainer::add_particle)
    .def("ene", [](qed::PhotonContainer& s) {return s.eneArr;}, py::return_value_policy::reference)
    // FIXME: use base class definition via getter/setter members to avoid this redefinition
    .def("loc",          [](qed::PhotonContainer& s, size_t idim) 
        {
          return s.loc(idim); 
        }, py::return_value_policy::reference)
    .def("vel",          [](qed::PhotonContainer& s, size_t idim) 
        {
          return s.vel(idim); 
        }, py::return_value_policy::reference);


  // 2D bindings
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  auto t2 = qed::declare_tile<2>(m_2d, "Tile");




}




} // ns qed
