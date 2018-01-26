#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../tools/mesh.h"
#include "amr/mesh.h"
#include "amr/numerics.h"
#include "amr/refiner.h"
#include "amr_momentum_solver.h"


typedef float Realf;
// typedef toolbox::AdaptiveMesh<Realf, 1> AM1d;
typedef toolbox::AdaptiveMesh<Realf, 3> AM3d;
typedef toolbox::Adapter<Realf, 3> Adapter3d;


//typedef std::vector<vlasov::PlasmaBlock> PlasmaBlockVector;
//PYBIND11_MAKE_OPAQUE(PlasmaBlockVector);




/// trampoline class for VlasovVelocitySolver
class PyMomentumSolver : public vlasov::MomentumSolver<Realf> {
  public:
    using vlasov::MomentumSolver<Realf>::MomentumSolver;
    using vlasov::MomentumSolver<Realf>::solve;

    void solveMesh( 
        AM3d& mesh0, 
        AM3d& mesh1, 
        std::array<Realf, 3>& E,
        std::array<Realf, 3>& B,
        Realf qm
        ) override {
      PYBIND11_OVERLOAD_PURE(
          void, 
          vlasov::MomentumSolver<Realf>, 
          solveMesh, 
          mesh0, mesh1, E, B, qm 
          );
    }
};





PYBIND11_MODULE(pyplasmaDev, m) {



  py::class_<AM3d >(m, "AdaptiveMesh3D")
    .def(py::init<>())
    .def_readwrite("length",                     &AM3d::length)
    .def_readwrite("maximum_refinement_level",   &AM3d::maximum_refinement_level)
    .def("resize",               &AM3d::resize)
    .def("get_cell_from_indices", &AM3d::get_cell_from_indices)
    .def("get_indices",          &AM3d::get_indices)
    .def("get_refinement_level", &AM3d::get_refinement_level)
    .def("get_parent_indices",   &AM3d::get_parent_indices)
    .def("get_parent",           &AM3d::get_parent)
    .def("get_maximum_possible_refinement_level",&AM3d::get_maximum_possible_refinement_level)
    .def("get_level_0_parent_indices", &AM3d::get_level_0_parent_indices)
    .def("get_level_0_parent",   &AM3d::get_level_0_parent)
    .def("get_children",         &AM3d::get_children)
    .def("get_siblings",         &AM3d::get_siblings)
    .def("get_cells",            &AM3d::get_cells)
    .def("__getitem__", [](const AM3d &s, py::tuple indx) 
        { 
        uint64_t i = indx[0].cast<uint64_t>();
        uint64_t j = indx[1].cast<uint64_t>();
        uint64_t k = indx[2].cast<uint64_t>();
        int    rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        return s.get_from_roots(cid);
        })
    .def("__setitem__", [](AM3d &s, py::tuple indx, Realf v) 
        { 
        uint64_t i = indx[0].cast<uint64_t>();
        uint64_t j = indx[1].cast<uint64_t>();
        uint64_t k = indx[2].cast<uint64_t>();
        int   rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        s.set(cid, v);
        })
    .def("clip_cells",       &AM3d::clip_cells)
    //.def("cut_roots",        &AM3d::cut_roots)
    // .def("split",                   &AM3d::split)
    // .def("update_from_children",    &AM3d::update_from_children)
    // .def("update_from_leafs",       &AM3d::update_from_leafs)
    // .def("get_from_leafs",          &AM3d::get_from_leafs)
    .def("is_leaf",                 &AM3d::is_leaf)
    .def("set_min",                 &AM3d::set_min)
    .def("set_max",                 &AM3d::set_max)
    .def("get_size",                &AM3d::get_size)
    .def("get_length",              &AM3d::get_length)
    .def("get_center",              &AM3d::get_center)
    .def("get_level_0_cell_length", &AM3d::get_level_0_cell_length);



    // -------------------------------------------------- 
    // AMR numerics
      
  m.def("deriv", &toolbox::deriv<Realf, 3>);

  m.def("grad",  &toolbox::grad<Realf, 3>);

  m.def("trilinear_interp", &toolbox::trilinear_interp<Realf>);

  m.def("tricubic_interp", &toolbox::tricubic_interp<Realf>);


  py::class_<Adapter3d>(m, "Adapter")
    .def(py::init<>())
    .def_readwrite("tolerance",         &Adapter3d::tolerance)
    .def_readwrite("cells_to_refine",   &Adapter3d::cells_to_refine)
    .def_readwrite("cells_to_unrefine", &Adapter3d::cells_to_unrefine)
    .def_readwrite("cells_created",     &Adapter3d::cells_created)
    .def_readwrite("cells_removed",     &Adapter3d::cells_removed)
    .def("set_maximum_data_value",      &Adapter3d::set_maximum_data_value)
    .def("maximum_value",               &Adapter3d::maximum_value)
    .def("maximum_gradient",            &Adapter3d::maximum_gradient)
    .def("check",                       &Adapter3d::check)
    .def("refine",                      &Adapter3d::refine)
    .def("unrefine",                    &Adapter3d::unrefine);


  py::class_<vlasov::MomentumSolver<Realf>, PyMomentumSolver> vvsol(m, "MomentumSolver");
  vvsol
    .def(py::init<>())
    .def("solve",     &vlasov::MomentumSolver<Realf>::solve)
    .def("solveMesh", &vlasov::MomentumSolver<Realf>::solveMesh);

  py::class_<vlasov::AmrMomentumLagrangianSolver<Realf>>(m, "AmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());




  py::class_<toolbox::Mesh<Realf,1>>(m, "Mesh")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &toolbox::Mesh<Realf,1>::Nx)
    .def_readwrite("Ny", &toolbox::Mesh<Realf,1>::Ny)
    .def_readwrite("Nz", &toolbox::Mesh<Realf,1>::Nz)
    .def("indx",         &toolbox::Mesh<Realf,1>::indx)
    .def("__getitem__", [](const toolbox::Mesh<Realf,1> &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < -1) throw py::index_error();
        if (j < -1) throw py::index_error();
        if (k < -1) throw py::index_error();

        if (i > (int)s.Nx+1) throw py::index_error();
        if (j > (int)s.Ny+1) throw py::index_error();
        if (k > (int)s.Nz+1) throw py::index_error();

        return s(i,j,k);
      })
    .def("__setitem__", [](toolbox::Mesh<Realf,1> &s, py::tuple indx, Realf val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < -1) throw py::index_error();
        if (j < -1) throw py::index_error();
        if (k < -1) throw py::index_error();

        if (i > (int)s.Nx+1) throw py::index_error();
        if (j > (int)s.Ny+1) throw py::index_error();
        if (k > (int)s.Nz+1) throw py::index_error();

        s(i,j,k) = val;
        })
    .def("clear",        &toolbox::Mesh<Realf,1>::clear);




  //typedef toolbox::Mesh<toolbox::AdaptiveMesh<Realf,3>,0> vmeshBlock;
  py::class_<vlasov::PlasmaBlock>(m, "PlasmaBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &vlasov::PlasmaBlock::Nx)
    .def_readwrite("Ny", &vlasov::PlasmaBlock::Ny)
    .def_readwrite("Nz", &vlasov::PlasmaBlock::Nz)
    .def_readwrite("qm", &vlasov::PlasmaBlock::qm)
    .def("__getitem__", [](const vlasov::PlasmaBlock &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        return s.block(i,j,k);
      }, py::return_value_policy::reference)
    .def("__setitem__", [](vlasov::PlasmaBlock &s, py::tuple indx, AM3d val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();
        int k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        s.block(i,j,k) = val;
        })
    .def("clear",       [](vlasov::PlasmaBlock &s){s.block.clear();});



};






