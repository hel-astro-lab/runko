#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#include "../coupling/tile.h"


namespace cpl {

//--------------------------------------------------
template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name) 
{

  return 
  py::class_<cpl::Tile<D>, 
             pic::Tile<D>,
             ffe::Tile<D>,
             emf::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<cpl::Tile<D>>
             >(m, 
               pyclass_name.c_str(),
               py::multiple_inheritance()
               )
    .def(py::init<int, int, int>())
    .def_readwrite("active_mode", &cpl::Tile<D>::active_mode);
}


// python bindings for coupled classes & functions
void bind_cpl(py::module& m_sub)
{

  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  auto t3 = cpl::declare_tile<3>(m_3d, "Tile_ffe_pic");

}

} // end of ns cpl
