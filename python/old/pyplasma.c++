#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../tools/mesh.h"

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


using Realf = float;
// typedef toolbox::AdaptiveMesh<Realf, 1> AM1d;
typedef toolbox::AdaptiveMesh<Realf, 3> AM3d;
typedef toolbox::Adapter<Realf, 3> Adapter3d;


/// trampoline class for VlasovVelocitySolver
typedef vlasov::MomentumSolver<Realf,3> momsol; // PYBIND preprocessor macro freaks out 
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
        Realf qm,
        Realf dt,
        Realf cfl
        ) override {
      PYBIND11_OVERLOAD_PURE(
          void, 
          momsol, 
          solve_mesh, 
          mesh0, mesh1, E, B, qm, dt, cfl
          );
    }
};


/// trampoline class for VlasovSpatialSolver
class PySpatialSolver : public vlasov::SpatialSolver<Realf> {
  public:

    void solve(
      vlasov::VlasovTile& tile,
      vlasov::Grid& grid
      ) override {
      PYBIND11_OVERLOAD_PURE(
          void,
          vlasov::SpatialSolver<Realf>,
          solve,
          tile, grid
          );
    }

};




  // generator for Mesh bindings with type T and halo H
  template<typename T, int H>
  void declare_mesh(
      py::module &m, 
      const std::string& pyclass_name) 
  {
      using Class = toolbox::Mesh<T, H>;
      py::class_<Class>(m, pyclass_name.c_str())

      .def(py::init<size_t, size_t, size_t>())
      .def_readwrite("Nx", &Class::Nx)
      .def_readwrite("Ny", &Class::Ny)
      .def_readwrite("Nz", &Class::Nz)
      .def("indx",         &Class::indx)
      .def("__getitem__", [](const Class &s, py::tuple indx) 
        {
          auto i = indx[0].cast<int>();
          auto j = indx[1].cast<int>();
          auto k = indx[2].cast<int>();


          // NOTE: these are out-of-bounds; not inbound checks
          if (i < -H) throw py::index_error();
          if (j < -H) throw py::index_error();
          if (k < -H) throw py::index_error();

          if (i >= (int)s.Nx+H) throw py::index_error();
          if (j >= (int)s.Ny+H) throw py::index_error();
          if (k >= (int)s.Nz+H) throw py::index_error();

          return s(i,j,k);
        })
      .def("__setitem__", [](Class &s, py::tuple indx, Realf val) 
        {
          auto i = indx[0].cast<int>();
          auto j = indx[1].cast<int>();
          auto k = indx[2].cast<int>();

          if (i < -H) throw py::index_error();
          if (j < -H) throw py::index_error();
          if (k < -H) throw py::index_error();

          if (i >= (int)s.Nx+H) throw py::index_error();
          if (j >= (int)s.Ny+H) throw py::index_error();
          if (k >= (int)s.Nz+H) throw py::index_error();

          s(i,j,k) = val;
          })
      .def("clear",        &Class::clear);
  }

  // generator for damped tile for various directions
  template<int S>
  void declare_PlasmaTileDamped(
      py::module& m,
      const std::string& pyclass_name) 
  {
    //using Class = fields::PlasmaTileDamped<S>; // does not function properly; maybe not triggering template?
    //have to use explicit name instead like this

    py::class_<fields::PlasmaTileDamped<S>,
             fields::PlasmaTile,
             std::shared_ptr<fields::PlasmaTileDamped<S>>
            >(m, pyclass_name.c_str() )
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t, size_t>())
    .def_readwrite("ex_ref",   &fields::PlasmaTileDamped<S>::ex_ref, py::return_value_policy::reference)
    .def_readwrite("ey_ref",   &fields::PlasmaTileDamped<S>::ey_ref, py::return_value_policy::reference)
    .def_readwrite("ez_ref",   &fields::PlasmaTileDamped<S>::ez_ref, py::return_value_policy::reference)
    .def_readwrite("bx_ref",   &fields::PlasmaTileDamped<S>::bx_ref, py::return_value_policy::reference)
    .def_readwrite("by_ref",   &fields::PlasmaTileDamped<S>::by_ref, py::return_value_policy::reference)
    .def_readwrite("bz_ref",   &fields::PlasmaTileDamped<S>::bz_ref, py::return_value_policy::reference)
    .def_readwrite("fld1",     &fields::PlasmaTileDamped<S>::fld1)
    .def_readwrite("fld2",     &fields::PlasmaTileDamped<S>::fld2)
    .def("damp_fields",         &fields::PlasmaTileDamped<S>::damp_fields);

  }




//--------------------------------------------------

// python bindings for plasma classes & functions
PYBIND11_MODULE(pyplasma, m) {

  // Loading tile bindings from corgi library
  py::object corgi_tile = (py::object) py::module::import("pycorgi").attr("Tile");


  /// General class for handling Maxwell's equations
  py::class_<fields::PlasmaTile,
             corgi::Tile, 
             std::shared_ptr<fields::PlasmaTile>
            >(m, "PlasmaTile")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",   &fields::PlasmaTile::dt)
    .def_readwrite("dx",   &fields::PlasmaTile::dx)
    .def_readwrite("cfl",  &fields::PlasmaTile::cfl)
    .def("cycle_yee",         &fields::PlasmaTile::cycle_yee)
    .def("cycle_current",     &fields::PlasmaTile::cycle_current)
    .def("cycle_current_2d",   &fields::PlasmaTile::cycle_current_2d)
    .def("push_e",            &fields::PlasmaTile::push_e)
    .def("push_half_b",        &fields::PlasmaTile::push_half_b)
    .def("deposit_current",   &fields::PlasmaTile::deposit_current)
    .def("get_yee",           &fields::PlasmaTile::get_yee, py::return_value_policy::reference)
    .def("get_analysis",      &fields::PlasmaTile::get_analysis, py::return_value_policy::reference)
    .def("add_analysis_species", &fields::PlasmaTile::add_analysis_species)
    .def("update_boundaries",  &fields::PlasmaTile::update_boundaries)
    .def("update_boundaries_2d",&fields::PlasmaTile::update_boundaries_2d)
    .def("exchange_currents",  &fields::PlasmaTile::exchange_currents)
    .def("exchange_currents_2d",&fields::PlasmaTile::exchange_currents_2d);



  declare_PlasmaTileDamped<-1>(m, "PlasmaTileDamped_LX");
  declare_PlasmaTileDamped<+1>(m, "PlasmaTileDamped_RX");
  declare_PlasmaTileDamped<-2>(m, "PlasmaTileDamped_LY");
  declare_PlasmaTileDamped<+2>(m, "PlasmaTileDamped_RY");


  // Loading node bindings from corgi library
  py::object corgi_node = (py::object) py::module::import("pycorgi").attr("Node");
  py::class_<vlasov::Grid>(m, "Grid", corgi_node)
    .def(py::init<size_t, size_t>());

  py::class_<fields::YeeLattice>(m, "YeeLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex",   &fields::YeeLattice::ex)
    .def_readwrite("ey",   &fields::YeeLattice::ey)
    .def_readwrite("ez",   &fields::YeeLattice::ez)
    .def_readwrite("bx",   &fields::YeeLattice::bx)
    .def_readwrite("by",   &fields::YeeLattice::by)
    .def_readwrite("bz",   &fields::YeeLattice::bz)
    .def_readwrite("jx",   &fields::YeeLattice::jx)
    .def_readwrite("jy",   &fields::YeeLattice::jy)
    .def_readwrite("jz",   &fields::YeeLattice::jz)
    .def_readwrite("jx1",  &fields::YeeLattice::jx1)
    .def_readwrite("rho",  &fields::YeeLattice::rho);


  py::class_<fields::PlasmaMomentLattice>(m, "PlasmaMomentLattice")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("rho",      &fields::PlasmaMomentLattice::rho)
    .def_readwrite("mgamma",   &fields::PlasmaMomentLattice::mgamma)
    .def_readwrite("Vx",       &fields::PlasmaMomentLattice::Vx)
    .def_readwrite("Vy",       &fields::PlasmaMomentLattice::Vy)
    .def_readwrite("Vz",       &fields::PlasmaMomentLattice::Vz)
    .def_readwrite("Tx",       &fields::PlasmaMomentLattice::Tx)
    .def_readwrite("Ty",       &fields::PlasmaMomentLattice::Ty)
    .def_readwrite("Tz",       &fields::PlasmaMomentLattice::Tz)
    .def_readwrite("ekin",     &fields::PlasmaMomentLattice::ekin);


  py::class_<vlasov::VlasovTile, 
             fields::PlasmaTile,
             corgi::Tile, 
             std::shared_ptr<vlasov::VlasovTile>
             >(m, "VlasovTile")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",     &vlasov::VlasovTile::dt)
    .def_readwrite("dx",     &vlasov::VlasovTile::dx)
    .def("get_plasma_species", [](vlasov::VlasovTile& tile, size_t i, size_t s) 
        { return tile.steps.get(i).at(s); }, py::return_value_policy::reference)
    .def("insert_initial_species", [](vlasov::VlasovTile& c, 
                                  std::vector<vlasov::PlasmaBlock> species){
        // push twice to initialize both time steps (current and future)
        c.steps.push_back(species);
        c.steps.push_back(species);

        c.Nspecies = species.size();

        })

    .def("clip",         &vlasov::VlasovTile::clip)
    .def("cycle",        &vlasov::VlasovTile::cycle);



  py::class_<AM3d >(m, "AdaptiveMesh3D")
    .def(py::init<>())
    .def_readwrite("length",                     &AM3d::length)
    .def_readwrite("maximum_refinement_level",   &AM3d::maximum_refinement_level)
    .def_readwrite("top_refinement_level",       &AM3d::top_refinement_level)
    .def("resize",               &AM3d::resize)
    .def("get_cell_from_indices", &AM3d::get_cell_from_indices)
    .def("get_indices",          &AM3d::get_indices)
    .def("get_refinement_level", &AM3d::get_refinement_level)
    .def("get_parent_indices",   &AM3d::get_parent_indices)
    .def("get_parent",           &AM3d::get_parent)
    .def("get_maximum_possible_refinement_level",&AM3d::get_maximum_possible_refinement_level)
    .def("set_maximum_refinement_level", &AM3d::set_maximum_refinement_level)
    .def("get_level_0_parent_indices", &AM3d::get_level_0_parent_indices)
    .def("get_level_0_parent",   &AM3d::get_level_0_parent)
    .def("get_children",         &AM3d::get_children)
    .def("get_siblings",         &AM3d::get_siblings)
    .def("get_cells",            &AM3d::get_cells)
    .def("__getitem__", [](const AM3d &s, py::tuple indx) 
        { 
        auto i = indx[0].cast<uint64_t>();
        auto j = indx[1].cast<uint64_t>();
        auto k = indx[2].cast<uint64_t>();
        auto    rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        return s.get_from_roots(cid);
        })
    .def("__setitem__", [](AM3d &s, py::tuple indx, Realf v) 
        { 
        auto i = indx[0].cast<uint64_t>();
        auto j = indx[1].cast<uint64_t>();
        auto k = indx[2].cast<uint64_t>();
        auto   rfl = indx[3].cast<int>();
        uint64_t cid = s.get_cell_from_indices({{i,j,k}}, rfl);

        if(cid == AM3d::error_cid) {throw py::index_error();}

        s.set(cid, v);
        })
    .def("clip_cells",              &AM3d::clip_cells)
    .def("is_leaf",                 &AM3d::is_leaf)
    .def("set_min",                 &AM3d::set_min)
    .def("set_max",                 &AM3d::set_max)
    .def("get_size",                &AM3d::get_size)
    .def("get_length",              &AM3d::get_length)
    .def("get_center",              &AM3d::get_center)
    .def("get_level_0_cell_length", &AM3d::get_level_0_cell_length);



    // -------------------------------------------------- 
    // AMR numerics
      
  //m.def("deriv", &toolbox::deriv<Realf, 3>);
  //m.def("grad",  &toolbox::grad<Realf, 3>);
  //m.def("interp_linear", &toolbox::interp_linear<Realf,3>);
  //m.def("interp_cubic", &toolbox::interp_cubic<Realf,3>);


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


  // general interface for momentum solvers
  py::class_<vlasov::MomentumSolver<Realf,3>, PyMomentumSolver > vvsol(m, "MomentumSolver");
  vvsol
    .def(py::init<>())
    .def("solve",     &vlasov::MomentumSolver<Realf,3>::solve)
    .def("solve_mesh", &vlasov::MomentumSolver<Realf,3>::solve_mesh);

  // AMR Lagrangian solver
  py::class_<vlasov::AmrMomentumLagrangianSolver<Realf,3>>(m, "AmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());

  py::class_<vlasov::GravityAmrMomentumLagrangianSolver<Realf,3>>(m, "GravityAmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());



  // general interface for spatial solvers
  py::class_<vlasov::SpatialSolver<Realf>, PySpatialSolver> vssol(m, "SpatialSolver");
  vssol
    .def(py::init<>())
    .def("solve", &vlasov::SpatialSolver<Realf>::solve);


  // AMR Lagrangian solver
  py::class_<vlasov::AmrSpatialLagrangianSolver<Realf>>(m, "AmrSpatialLagrangianSolver", vssol)
    .def(py::init<>());


  /// Vlasov tile analyzator
  py::class_<vlasov::Analyzator<Realf> >(m, "Analyzator")
    .def(py::init<>());


  declare_mesh<Realf, 0>(m, std::string("Mesh0") );
  declare_mesh<Realf, 1>(m, std::string("Mesh1") );
  declare_mesh<Realf, 3>(m, std::string("Mesh3") );



  py::class_<vlasov::PlasmaBlock>(m, "PlasmaBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("Nx", &vlasov::PlasmaBlock::Nx)
    .def_readwrite("Ny", &vlasov::PlasmaBlock::Ny)
    .def_readwrite("Nz", &vlasov::PlasmaBlock::Nz)
    .def_readwrite("qm", &vlasov::PlasmaBlock::qm)
    .def("__getitem__", [](const vlasov::PlasmaBlock &s, py::tuple indx) 
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
    .def("__setitem__", [](vlasov::PlasmaBlock &s, py::tuple indx, AM3d val) 
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
    .def("clear",       [](vlasov::PlasmaBlock &s){s.block.clear();});


    m.def("initial_step_1d", &vlasov::initial_step<1>);
    m.def("initial_step", &vlasov::initial_step<3>);

    m.def("step_location", &vlasov::step_location);


    m.def("step_velocity_1d", &vlasov::step_velocity<1>);
    m.def("step_velocity",   &vlasov::step_velocity<3>);
    m.def("step_velocity_with_gravity_1d",   &vlasov::step_velocity_with_gravity<1>);

    m.def("analyze",      &vlasov::analyze);

    m.def("write_yee", &vlasov::write_yee);

    m.def("write_analysis", &vlasov::write_analysis);

    m.def("write_mesh", &vlasov::write_mesh);

}



