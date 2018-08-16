#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../pic/boundaries/wall.h"



// generator for wall cell
template<int S>
void declare_PicCellWall(
    py::module& m,
    std::string pyclass_name) 
{

  py::class_<pic::PicCellWall<S>,
    fields::PlasmaCellDamped<S>,
    pic::PicCell,
    std::shared_ptr<pic::PicCellWall<S>>
      >(m, pyclass_name.c_str() )
      .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>());

}




PYBIND11_MODULE(pypic_boundaries, m) {

  py::object corgiCell  = (py::object) py::module::import("pycorgi").attr("Cell");
  py::object vlasovCell = (py::object) py::module::import("pyplasma").attr("VlasovCell");


  declare_PicCellWall<-1>(m, "PicCellWall_LX");
  declare_PicCellWall<+1>(m, "PicCellWall_RX");

  declare_PicCellWall<-2>(m, "PicCellWall_LY");
  declare_PicCellWall<+2>(m, "PicCellWall_RY");

}









