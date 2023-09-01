#include "py_submodules.h"

#include <pybind11/numpy.h>

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

  tuple_pair get_minmax_ene( string t1, string t2, double ene) override { 
    PYBIND11_OVERLOAD_PURE(
        tuple_pair, // return type
        Interaction,                // parent class
        get_minmax_ene,             // name of function in C++
        t1, // arguments
        t2,
        ene
        );
  }

  tuple_pair comp_cross_section( 
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2) override {
    PYBIND11_OVERLOAD_PURE(
        tuple_pair, // return type
        Interaction,                  // parent class
        comp_cross_section,           // name of function in C++
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

  tuple_pair accumulate( string t1, float_p e1, string t2, float_p e2) override { 
    PYBIND11_OVERLOAD_PURE(
        tuple_pair, // return type
        Interaction,               // parent class
        accumulate,                // name of function in C++
        t1, // arguments
        e1,
        t2,
        e2
        );
  }

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
    .def_readwrite("do_accumulate", &qed::Interaction::do_accumulate)
    .def(py::init<string, string>())
    .def("get_minmax_ene",      &qed::Interaction::get_minmax_ene)
    .def("comp_cross_section",  &qed::Interaction::comp_cross_section)
    .def("accumulate",          &qed::Interaction::accumulate)
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

  // Compton  
  py::class_<qed::Compton>(m_sub, "Compton", qedinter)
    .def_readwrite("no_electron_update", &qed::Compton::no_electron_update)
    .def_readwrite("no_photon_update",   &qed::Compton::no_photon_update)
    .def_readwrite("ming",               &qed::Compton::ming)
    .def_readwrite("minx2z",             &qed::Compton::minx2z)
    .def(py::init<string, string>());



  py::module m_2d = m_sub.def_submodule("twoD",   "2D specializations");
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  //--------------------------------------------------
  // Particle pairing routines
  py::class_<qed::Pairing<2> >(m_2d, "Pairing")
    .def(py::init<>())
    .def_readwrite("prob_norm",&qed::Pairing<2>::prob_norm)
    .def("solve",              &qed::Pairing<2>::solve)
    .def("solve_mc",           &qed::Pairing<2>::solve_mc)
    .def("add_interaction",    &qed::Pairing<2>::add_interaction, py::keep_alive<1,2>() )
    .def("rescale",            &qed::Pairing<2>::rescale);


  py::class_<qed::Pairing<3> >(m_3d, "Pairing")
    .def(py::init<>())
    .def_readwrite("prob_norm",   &qed::Pairing<3>::prob_norm)
    .def_readwrite("leaked_ene",  &qed::Pairing<3>::leaked_ene)
    .def_readwrite("leaked_wsum", &qed::Pairing<3>::leaked_wsum)
    .def_readwrite("leaked_pnum", &qed::Pairing<3>::leaked_pnum)
    .def_readwrite("inj_ene_ph",  &qed::Pairing<3>::inj_ene_ph)
    .def_readwrite("inj_ene_ep",  &qed::Pairing<3>::inj_ene_ep)
    .def_readwrite("tau_measured",&qed::Pairing<3>::tau_measured)
    .def("add_interaction",    &qed::Pairing<3>::add_interaction, py::keep_alive<1,2>() )
    .def("rescale",            &qed::Pairing<3>::rescale)
    .def("inject_photons",     &qed::Pairing<3>::inject_photons)
    .def("inject_plaw_pairs",  &qed::Pairing<3>::inject_plaw_pairs)
    .def("leak_photons",       &qed::Pairing<3>::leak_photons)
    .def("update_hist_lims",   &qed::Pairing<3>::update_hist_lims)
    .def("clear_hist",         &qed::Pairing<3>::clear_hist)
    .def("solve",              &qed::Pairing<3>::solve)
    .def("solve_mc",           &qed::Pairing<3>::solve_mc)
    .def("get_hist_edges",   [](qed::Pairing<3>& s)
        {
          const auto N = static_cast<pybind11::ssize_t>(s.hist_nbin);
          auto v = pybind11::array_t<double>( {N}, s.hist_ene_edges.data() );
          return v;
        })
    .def("get_hist_cnts",    [](qed::Pairing<3>& s)
        {
          const auto N = static_cast<pybind11::ssize_t>(s.hist_nbin);
          auto v = pybind11::array_t<double>( {N}, s.hist.data() );
          return v;
        });

}




} // ns qed
