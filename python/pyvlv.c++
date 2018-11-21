#include "py_submodules.h"


#include "../definitions.h"
//#include "../tools/mesh.h"


//--------------------------------------------------
// Vlasov module
  
#include "../vlasov/tile.h"
//#include "../vlasov/grid.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

#include "../vlasov/amr/mesh.h"
//#include "../vlasov/amr/numerics.h"
#include "../vlasov/amr/refiner.h"
#include "../vlasov/amr/operators.h"

#include "../vlasov/momentum-solvers/amr_momentum_solver.h"
#include "../vlasov/momentum-solvers/bwd_lagrangian.h"
#include "../vlasov/momentum-solvers/bwd_lagrangian_gravity.h"

#include "../vlasov/spatial-solvers/amr_spatial_solver.h"
#include "../vlasov/amr_analyzator.h"

#include "../vlasov/tasker.h"


// boundaries
#include "../vlasov/boundaries/outflow.h"
#include "../vlasov/boundaries/piston.h"



namespace vlv {


using Adapter3d = toolbox::Adapter<Realf, 3>;
using AM1d = toolbox::AdaptiveMesh<Realf, 1>;
using AM3d = toolbox::AdaptiveMesh<Realf, 3>;

//--------------------------------------------------
/// trampoline class for VlasovVelocitySolver
using momsol = vlv::MomentumSolver<Realf,1,1>; // PYBIND preprocessor macro freaks out 
                                             // of commas so we hide them with typedef

                                                  
class PyMomentumSolver : public momsol {

  public:
    using momsol::MomentumSolver;
    using momsol::solve;

    void solve_mesh( 
        AM3d& mesh0, 
        AM3d& mesh1, 
        std::array<Realf, 3>& E,
        std::array<Realf, 3>& B,
        vlv::tools::Params<Realf>& params
        ) override {
      PYBIND11_OVERLOAD_PURE(
          void, 
          momsol, 
          solve_mesh, 
          mesh0, mesh1, E, B, params 
          );
    }
};


/// trampoline class for VlasovSpatialSolver
class PySpatialSolver : public vlv::SpatialSolver<Realf> {
  public:

    void solve(
      vlv::Tile<1>& tile,
      corgi::Node<1>& grid
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

//template<size_t D>
//auto declare_grid(
//    py::module& m,
//    const std::string& pyclass_name,
//    const std::string& pycorgi_name
//    ) 
//{
//  py::object corgi_node = 
//    (py::object) py::module::import("pycorgi").attr(pycorgi_name.c_str());
//
//  return 
//    py::class_<vlv::Grid<D> >(m, pyclass_name.c_str(), corgi_node);
//      //.def()
//}


//--------------------------------------------------
// generator for vlv::Tile 

template<size_t D>
auto declare_tile(
    py::module& m,
    const std::string& pyclass_name
    ) 
{

  return
    py::class_<vlv::Tile<D>, 
             fields::Tile<D>,
             std::shared_ptr<vlv::Tile<D> >
             >(m, pyclass_name.c_str())
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("dx",        &vlv::Tile<D>::dx)
    .def_readwrite("threshold", &vlv::Tile<D>::threshold)
    .def("get_plasma_species", [](vlv::Tile<D>& tile, size_t i, size_t s) 
        { return tile.steps.get(i).at(s); }, py::return_value_policy::reference)
    .def("insert_initial_species", [](vlv::Tile<D>& c, 
                                    std::vector<vlv::PlasmaBlock> species){
        // push twice to initialize both time steps (current and future)
        c.steps.push_back(species);
        c.steps.push_back(species);

        c.Nspecies = species.size();

        })

    .def("clip",           &vlv::Tile<D>::clip)
    .def("clip_neighbors", &vlv::Tile<D>::clip_neighbors)
    .def("cycle",          &vlv::Tile<D>::cycle);

}


//--------------------------------------------------
namespace outflow {
  // generator for outflow tile
  template<size_t D, int S>
    auto declare_tile(
        py::module& m,
        const std::string& pyclass_name) 
    {
      return
        py::class_<
        vlv::outflow::Tile<D, S>,
          fields::damping::Tile<D, S>,
          vlv::Tile<D>,
        std::shared_ptr<vlv::outflow::Tile<D,S>>
        >(m, pyclass_name.c_str() )
    .def(py::init<size_t, size_t, size_t>());
    }
}


namespace piston {
  // generator for piston tile
  template<size_t D>
    auto declare_tile(
        py::module& m,
        const std::string& pyclass_name) 
    {
      return
        py::class_<vlv::piston::Tile<D>,
        vlv::Tile<D>,
      std::shared_ptr<vlv::piston::Tile<D>>
        >(m, pyclass_name.c_str() )
    .def(py::init<size_t, size_t, size_t>())
    .def("reflect", &vlv::piston::Tile<D>::reflect);
    }
}



/// python bindings for Vlasov module
void bind_vlv(py::module& m_sub)
{
  //--------------------------------------------------
  // general bindings

  py::class_<vlv::PlasmaBlock>(m_sub, "PlasmaBlock")
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
  // 1D bindings

  py::module m_1d = m_sub.def_submodule("oneD", "1D specializations");

  //--------------------------------------------------
  // Tile bindings
    
  auto t1 = declare_tile<1>(m_1d, "Tile");


  //--------------------------------------------------

  py::class_<Adapter3d>(m_1d, "Adapter")
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
  py::class_<vlv::MomentumSolver<Realf,1,1>, PyMomentumSolver > vvsol(m_1d, "MomentumSolver");
  vvsol
    .def(py::init<>())
    .def("solve",     &vlv::MomentumSolver<Realf,1,1>::solve)
    .def("solve_mesh", &vlv::MomentumSolver<Realf,1,1>::solve_mesh);

  // AMR Lagrangian solver
  py::class_<vlv::AmrMomentumLagrangianSolver<Realf,1,1>>(m_1d, "AmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());

  py::class_<vlv::GravityAmrMomentumLagrangianSolver<Realf,1,1>>(m_1d, "GravityAmrMomentumLagrangianSolver", vvsol)
     .def(py::init<Realf,Realf>());


  //--------------------------------------------------

  // general interface for spatial solvers
  py::class_<vlv::SpatialSolver<Realf>, PySpatialSolver> vssol(m_1d, "SpatialSolver");
  vssol
    .def(py::init<>())
    .def("solve", &vlv::SpatialSolver<Realf>::solve);


  // AMR Lagrangian solver
  py::class_<vlv::AmrSpatialLagrangianSolver<Realf>>(m_1d, "AmrSpatialLagrangianSolver", vssol)
    .def(py::init<>());


  //--------------------------------------------------

  /// Vlasov tile analyzator
  py::class_<vlv::Analyzator<Realf> >(m_1d, "Analyzator")
    .def(py::init<>());


  //--------------------------------------------------


  m_1d.def("initial_step",    &vlv::initial_step<1>);
  m_1d.def("step_location",   &vlv::step_location);
  m_1d.def("step_velocity",   &vlv::step_velocity<1>);
  m_1d.def("step_velocity_with_gravity",   &vlv::step_velocity_with_gravity<1>);

  m_1d.def("analyze",        &vlv::analyze);
  m_1d.def("write_yee",       &vlv::write_yee<1>);
  m_1d.def("write_analysis",  &vlv::write_analysis<1>);
  m_1d.def("write_mesh",      &vlv::write_mesh<1>);

  m_1d.def("read_yee",        &vlv::read_yee<1>);
  m_1d.def("read_mesh",       &vlv::read_mesh<1>);


  //--------------------------------------------------
  // boundaries
  auto tf_R = outflow::declare_tile<1,+1>(m_1d, "Tile_outflow_L");
  auto tf_L = outflow::declare_tile<1,-1>(m_1d, "Tile_outflow_R");

  //auto tp = piston::declare_tile<1>(m_1d, "Tile_piston");
  



  //--------------------------------------------------
  // 2D bindings

  py::module m_2d = m_sub.def_submodule("twoD", "2D specializations");
  m_2d.def("write_yee",       &vlv::write_yee<2>);
  m_2d.def("write_analysis",  &vlv::write_analysis<2>);

  m_2d.def("read_yee",        &vlv::read_yee<2>);


}

} // end of namespace vlv
