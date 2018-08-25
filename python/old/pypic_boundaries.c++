#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../pic/boundaries/wall.h"



// generator for wall tile
template<int S>
void declare_PicTileWall(
    py::module& m,
    const std::string& pyclass_name) 
{

  py::class_<pic::PicTileWall<S>,
    fields::PlasmaTileDamped<S>,
    pic::PicTile,
    std::shared_ptr<pic::PicTileWall<S>>
      >(m, pyclass_name.c_str() )
      .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>());

}




PYBIND11_MODULE(pypic_boundaries, m) {

  py::object corgiTile  = (py::object) py::module::import("pycorgi").attr("Tile");
  py::object vlasovTile = (py::object) py::module::import("pyplasma").attr("VlasovTile");


  declare_PicTileWall<-1>(m, "PicTileWall_LX");
  declare_PicTileWall<+1>(m, "PicTileWall_RX");

  declare_PicTileWall<-2>(m, "PicTileWall_LY");
  declare_PicTileWall<+2>(m, "PicTileWall_RY");

}









