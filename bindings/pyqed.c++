#include "py_submodules.h"

//#include "../qed/photon.h"
//#include "../qed/tile.h"

#include "../qed/interactions/interaction.h"
#include "../qed/interactions/pair_ann.h"
#include "../qed/interactions/phot_ann.h"
#include "../qed/interactions/compton.h"
#include "../qed/pairing.h"



namespace qed {

  using std::string;

//--------------------------------------------------
//template<size_t D>
//auto declare_tile(
//    py::module& m,
//    const std::string& pyclass_name) 
//{
//
//  return 
//  py::class_<qed::Tile<D>,
//             pic::Tile<D>, 
//             fields::Tile<D>,
//             corgi::Tile<D>, 
//             std::shared_ptr<qed::Tile<D>>
//             >(m, 
//               pyclass_name.c_str(),
//               py::multiple_inheritance()
//               )
//    .def(py::init<int, int, int>())
//    .def("get_bucket",       &qed::Tile<D>::get_bucket, 
//        py::return_value_policy::reference)
//    .def("push_back",       &qed::Tile<D>::push_back);
//}

/// trampoline class for Interactions 
class PyInteraction : public Interaction
{
public:


  using Interaction::Interaction;


  // trampoline for each  virtual function
 
  using tuple_pair = std::tuple<float_p, float_p>; // following macro does not accept commas so we define this

  tuple_pair get_minmax_ene( string t1, string t2 ) override { 
    PYBIND11_OVERLOAD_PURE(
        tuple_pair, // return type
        Interaction,                // parent class
        get_minmax_ene,             // name of function in C++
        t1, // arguments
        t2  
        );
  }

  float_p comp_cross_section( 
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2) override {
    PYBIND11_OVERLOAD_PURE(
        float_p,                     // return type
        Interaction,                // parent class
        comp_cross_section,         // name of function in C++
        t1, ux1, uy1, uz1, t2, ux2, uy2, uz2
        );
  }

  //std::tuple<
  //  string, float_p, float_p, float_p,
  //  string, float_p, float_p, float_p>
  //    interact(
  //      string t1, float_p ux1, float_p uy1, float_p uz1,
  //      string t2, float_p ux2, float_p uy2, float_p uz2) override {
  //  PYBIND11_OVERLOAD_PURE(
  //      std::tuple<string, float_p, float_p, float_p, string, float_p, float_p, float_p>, // return type
  //      Interaction,                // parent class
  //      interact,                   // name of function in C++
  //      string, // arguments
  //      float_p, float_p, float_p,
  //      string,
  //      float_p, float_p, float_p
  //      );
  //  }

  void interact(
        string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
        string& t2, float_p& ux2, float_p& uy2, float_p& uz2) override {
    PYBIND11_OVERLOAD_PURE(
        void,                       // return type
        Interaction,                // parent class
        interact,                   // name of function in C++
        t1, ux1, uy1, uz1, t2, ux2, uy2, uz2
        );
    }

};


//--------------------------------------------------

// python bindings for radiation classes & functions
void bind_qed(py::module& m_sub)
{

  //py::class_<qed::PhotonContainer, pic::ParticleContainer<3>>(m_sub, "PhotonContainer")
  //  .def(py::init<>())
  //  .def("add_particle",  (void (qed::PhotonContainer::*)
  //        ( std::vector<float_p>, std::vector<float_p>, float_p, float_p ) ) 
  //          &qed::PhotonContainer::add_particle)
  //  .def("ene", [](qed::PhotonContainer& s) {return s.eneArr;}, py::return_value_policy::reference)
  //  // FIXME: use base class definition via getter/setter members to avoid this redefinition
  //  .def("loc",          [](qed::PhotonContainer& s, size_t idim) 
  //      {
  //        return s.loc(idim); 
  //      }, py::return_value_policy::reference)
  //  .def("vel",          [](qed::PhotonContainer& s, size_t idim) 
  //      {
  //        return s.vel(idim); 
  //      }, py::return_value_policy::reference);


  // 2D bindings
  //py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  //auto t2 = qed::declare_tile<2>(m_2d, "Tile");



  py::class_< qed::Interaction, std::shared_ptr<qed::Interaction>, PyInteraction > qedinter(m_sub, "Interaction");
  qedinter
    //.def(py::init<string&, string& >())
    .def(py::init<string, string>())
    //.def("get_ene",      &qed::Interaction::get_ene);
    .def("get_minmax_ene",      &qed::Interaction::get_minmax_ene)
    .def("comp_cross_section",  &qed::Interaction::comp_cross_section)
    //.def("interact",            &qed::Interaction::interact);
    .def("interact", [](qed::Interaction &self, 
          std::string t1, float_p ux1, float_p uy1, float_p uz1,
          std::string t2, float_p ux2, float_p uy2, float_p uz2) 
        {
          self.interact(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2); 
          return std::make_tuple(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2);
        });

  // Pair annihilation 
  py::class_<qed::PairAnn>(m_sub, "PairAnn", qedinter)
    .def(py::init<string, string>());

  // Photon annihilation 
  py::class_<qed::PhotAnn>(m_sub, "PhotAnn", qedinter)
    .def(py::init<string, string>());

  // Photon annihilation 
  py::class_<qed::Compton>(m_sub, "Compton", qedinter)
    .def(py::init<string, string>());


  //--------------------------------------------------
  // Particle pairing routines
  py::class_<qed::Pairing>(m_sub, "Pairing")
    .def(py::init<>())
    .def("solve",           &qed::Pairing::solve<3>)
    .def("add_interaction", &qed::Pairing::add_interaction, py::keep_alive<1,2>() );



}




} // ns qed
