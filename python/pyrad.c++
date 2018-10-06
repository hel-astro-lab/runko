#include "py_submodules.h"

#include "../radiation/photon.h"
#include "../radiation/tile.h"


namespace rad {

//--------------------------------------------------
template<size_t D>
auto declare_Tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<rad::Tile<D>,
             pic::Tile<D>, 
             fields::Tile<D>,
             std::shared_ptr<rad::Tile<D>>
             >(m, pyclass_name.c_str())
    .def(py::init<size_t, size_t, size_t>())
    .def("get_bucket",       &rad::Tile<D>::get_bucket, 
        py::return_value_policy::reference)
    .def("push_back",       &rad::Tile<D>::push_back);
}



//--------------------------------------------------

// python bindings for radiation classes & functions
void bind_rad(py::module& m_sub)
{

  // example of binding polymorphic classes
  //py::class_<Pet>(m, "Pet")
  //   .def(py::init<const std::string &, int>())
  //   .def("set", (void (Pet::*)(int)) &Pet::set, "Set the pet's age")
  //   .def("set", (void (Pet::*)(const std::string &)) &Pet::set, "Set the pet's name");

  py::class_<rad::PhotonBlock, pic::ParticleBlock>(m_sub, "PhotonBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def("add_particle",  (void (rad::PhotonBlock::*)
          ( std::vector<Realf>, std::vector<Realf>, Realf, Realf ) ) 
            &rad::PhotonBlock::add_particle)
    .def("wgt", [](rad::PhotonBlock& s) {return s.wgtArr;}, py::return_value_policy::reference)
    .def("ene", [](rad::PhotonBlock& s) {return s.eneArr;}, py::return_value_policy::reference)
    // FIXME: use base class definition via getter/setter members to avoid this redefinition
    .def("loc",          [](rad::PhotonBlock& s, size_t idim) 
        {
          return s.loc(idim); 
        }, py::return_value_policy::reference)
    .def("vel",          [](rad::PhotonBlock& s, size_t idim) 
        {
          return s.vel(idim); 
        }, py::return_value_policy::reference);


  // 2D bindings
  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  auto t2 = rad::declare_Tile<2>(m_2d, "Tile");




};




} // ns rad
