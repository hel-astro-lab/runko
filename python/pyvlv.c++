#include "py_submodules.h"


#include "../definitions.h"
//#include "../tools/mesh.h"


//--------------------------------------------------
// Vlasov module
  
//#include "../vlasov/tile.h"
#include "../vlasov/grid.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

#include "../vlasov/amr/mesh.h"
//#include "../vlasov/amr/numerics.h"
#include "../vlasov/amr/refiner.h"
#include "../vlasov/amr/operators.h"
#include "../vlasov/amr_momentum_solver.h"
#include "../vlasov/amr_spatial_solver.h"
#include "../vlasov/amr_analyzator.h"

#include "../vlasov/tasker.h"



using Adapter3d = toolbox::Adapter<Realf, 3>;
using AM1d = toolbox::AdaptiveMesh<Realf, 1>;
using AM3d = toolbox::AdaptiveMesh<Realf, 3>;

//--------------------------------------------------
/// trampoline class for VlasovVelocitySolver
using momsol = vlv::MomentumSolver<Realf,3>; // PYBIND preprocessor macro freaks out 
                                             // of commas so we hide them with typedef

                                                  
class PyMomentumSolver : public momsol {

  public:
    using momsol::MomentumSolver;
    using momsol::solve;

    void solveMesh( 
        AM3d& mesh0, 
        AM3d& mesh1, 
        std::array<Realf, 3>& E,
        std::array<Realf, 3>& B,
        Realf qm,
        Realf cfl
        ) override {
      PYBIND11_OVERLOAD_PURE(
          void, 
          momsol, 
          solveMesh, 
          mesh0, mesh1, E, B, qm, cfl
          );
    }
};


/// trampoline class for VlasovSpatialSolver
class PySpatialSolver : public vlv::SpatialSolver<Realf> {
  public:

    void solve(
      vlv::Tile<1>& tile,
      vlv::Grid<1>& grid
      ) override {
      PYBIND11_OVERLOAD_PURE(
          void,
          vlv::SpatialSolver<Realf>,
          solve,
          tile, grid
          );
    }
};


//--------------------------------------------------
// generator for vlv::Grid 

template<size_t D>
auto declare_Grid(
    py::module& m,
    const std::string& pyclass_name,
    const std::string& pycorgi_name
    ) 
{
  py::object corgiNode = 
    (py::object) py::module::import("pycorgi").attr(pycorgi_name.c_str());

  return 
    py::class_<vlv::Grid<D> >(m, pyclass_name.c_str(), corgiNode);
      //.def()
}


//--------------------------------------------------
// generator for vlv::Tile 

template<size_t D>
auto declare_Tile(
    py::module& m,
    const std::string& pyclass_name
    ) 
{

  return
    py::class_<vlv::Tile<D>, 
             fields::Tile<D>,
             corgi::Tile<D>, 
             std::shared_ptr<vlv::Tile<D> >
             >(m, pyclass_name.c_str())
    .def_readwrite("dx",     &vlv::Tile<D>::dx)
    .def("getPlasmaSpecies", [](vlv::Tile<D>& tile, size_t i, size_t s) 
        { return tile.steps.get(i).at(s); }, py::return_value_policy::reference)
    .def("insertInitialSpecies", [](vlv::Tile<D>& c, 
                                    std::vector<vlv::PlasmaBlock> species){
        // push twice to initialize both time steps (current and future)
        c.steps.push_back(species);
        c.steps.push_back(species);

        c.Nspecies = species.size();

        })

    .def("clip",         &vlv::Tile<D>::clip)
    .def("cycle",        &vlv::Tile<D>::cycle);

}







void bind_vlv(py::module& m)
{

  //--------------------------------------------------
  // Grid bindings

  //auto g1 = declare_Grid<1>(m, "Grid1D", "Node1D");
  auto g2 = declare_Grid<2>(m, "Grid2D", "Node");
  //auto g3 = declare_Grid<3>(m, "Grid3D", "Node3D");

  //g1.def(py::init<size_t>());
  g2.def(py::init<size_t, size_t>());
  //g3.def(py::init<size_t, size_t, size_t>());


  //--------------------------------------------------
  // Tile bindings
    
  auto t1 = declare_Tile<1>(m, "Tile1D");
  //auto t2 = declare_Tile<2>(m, "Tile2D");
  ////auto t3 = declare_Tile<3>(m, "Tile3D");

  t1.def(py::init<size_t>());
  //t2.def(py::init<size_t, size_t>());
  //t3.def(py::init<size_t, size_t, size_t>());




  //--------------------------------------------------

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

  //--------------------------------------------------

  // general interface for momentum solvers
  py::class_<vlv::MomentumSolver<Realf,3>, PyMomentumSolver > vvsol(m, "MomentumSolver");
  vvsol
    .def(py::init<>())
    .def("solve",     &vlv::MomentumSolver<Realf,3>::solve)
    .def("solveMesh", &vlv::MomentumSolver<Realf,3>::solveMesh);

  // AMR Lagrangian solver
  py::class_<vlv::AmrMomentumLagrangianSolver<Realf,3>>(m, "AmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());

  py::class_<vlv::GravityAmrMomentumLagrangianSolver<Realf,3>>(m, "GravityAmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());


  //--------------------------------------------------

  // general interface for spatial solvers
  py::class_<vlv::SpatialSolver<Realf>, PySpatialSolver> vssol(m, "SpatialSolver");
  vssol
    .def(py::init<>())
    .def("solve", &vlv::SpatialSolver<Realf>::solve);


  // AMR Lagrangian solver
  py::class_<vlv::AmrSpatialLagrangianSolver<Realf>>(m, "AmrSpatialLagrangianSolver", vssol)
    .def(py::init<>());


  //--------------------------------------------------

  /// Vlasov tile analyzator
  py::class_<vlv::Analyzator<Realf> >(m, "Analyzator")
    .def(py::init<>());


  //--------------------------------------------------

  py::class_<vlv::PlasmaBlock>(m, "PlasmaBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &vlv::PlasmaBlock::Nx)
    .def_readwrite("Ny", &vlv::PlasmaBlock::Ny)
    .def_readwrite("Nz", &vlv::PlasmaBlock::Nz)
    .def_readwrite("qm", &vlv::PlasmaBlock::qm)
    .def("__getitem__", [](const vlv::PlasmaBlock &s, py::tuple indx) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        return s.block(i,j,k);
      }, py::return_value_policy::reference)
    .def("__setitem__", [](vlv::PlasmaBlock &s, py::tuple indx, AM3d val) 
      {
        auto i = indx[0].cast<int>();
        auto j = indx[1].cast<int>();
        auto k = indx[2].cast<int>();

        if (i < 0) throw py::index_error();
        if (j < 0) throw py::index_error();
        if (k < 0) throw py::index_error();

        if (i > (int)s.Nx) throw py::index_error();
        if (j > (int)s.Ny) throw py::index_error();
        if (k > (int)s.Nz) throw py::index_error();

        s.block(i,j,k) = val;
        })
    .def("clear",       [](vlv::PlasmaBlock &s){s.block.clear();});



  //--------------------------------------------------


  m.def("stepInitial1d",  &vlv::stepInitial<1>);
  //m.def("stepInitial",    &vlv::stepInitial<3>);

  m.def("stepLocation",   &vlv::stepLocation);


  m.def("stepVelocity1d", &vlv::stepVelocity<1>);
  //m.def("stepVelocity",   &vlv::stepVelocity<3>);
  m.def("stepVelocityGravity1d",   &vlv::stepVelocityGravity<1>);

  m.def("analyze",        &vlv::analyze);

  m.def("writeYee",       &vlv::writeYee);

  m.def("writeAnalysis",  &vlv::writeAnalysis);

  m.def("writeMesh",      &vlv::writeMesh);






}
