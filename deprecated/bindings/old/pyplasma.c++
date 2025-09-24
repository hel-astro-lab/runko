#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include "../definitions.h"
#include "../tools/mesh.h"

//--------------------------------------------------
// Vlasov module
  
//#include "../vlv/tile.h"
#include "../vlv/grid.h"

#include "../em-fields/tile.h"
#include "../em-fields/damping_tile.h"

#include "../vlv/amr/mesh.h"
//#include "../vlv/amr/numerics.h"
#include "../vlv/amr/refiner.h"
#include "../vlv/amr/operators.h"
#include "../vlv/amr_momentum_solver.h"
#include "../vlv/amr_spatial_solver.h"
#include "../vlv/amr_analyzator.h"

#include "../vlv/tasker.h"


using float_m = float;
// typedef toolbox::AdaptiveMesh<float_m, 1> AM1d;
typedef toolbox::AdaptiveMesh<float_m, 3> AM3d;
typedef toolbox::Adapter<float_m, 3> Adapter3d;


/// trampoline class for VlasovVelocitySolver
typedef vlv::MomentumSolver<float_m,3> momsol; // PYBIND preprocessor macro freaks out 
                                                // of commas so we hide them with typedef
                                                  
class PyMomentumSolver : public momsol {

  public:
    using momsol::MomentumSolver;
    using momsol::solve;

    void solve_mesh( 
        AM3d& mesh0, 
        AM3d& mesh1, 
        std::array<float_m, 3>& E,
        std::array<float_m, 3>& B,
        float_m qm,
        float_m dt,
        float_m cfl
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
class PySpatialSolver : public vlv::SpatialSolver<float_m> {
  public:

    void solve(
      vlv::VlasovTile& tile,
      vlv::Grid& grid
      ) override {
      PYBIND11_OVERLOAD_PURE(
          void,
          vlv::SpatialSolver<float_m>,
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
      .def("__setitem__", [](Class &s, py::tuple indx, float_m val) 
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
    .def("cycle_grids",         &fields::PlasmaTile::cycle_grids)
    .def("cycle_current",     &fields::PlasmaTile::cycle_current)
    .def("cycle_current_2d",   &fields::PlasmaTile::cycle_current_2d)
    .def("push_e",            &fields::PlasmaTile::push_e)
    .def("push_half_b",        &fields::PlasmaTile::push_half_b)
    .def("deposit_current",   &fields::PlasmaTile::deposit_current)
    .def("get_grids",           &fields::PlasmaTile::get_grids, py::return_value_policy::reference)
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


  // Loading grid bindings from corgi library
  py::object corgi_node = (py::object) py::module::import("pycorgi").attr("Grid");
  py::class_<vlv::Grid>(m, "Grid", corgi_node)
    .def(py::init<size_t, size_t>());

  py::class_<fields::Grids>(m, "Grids")
    .def(py::init<size_t, size_t, size_t>())
    .def_readwrite("ex",   &fields::Grids::ex)
    .def_readwrite("ey",   &fields::Grids::ey)
    .def_readwrite("ez",   &fields::Grids::ez)
    .def_readwrite("bx",   &fields::Grids::bx)
    .def_readwrite("by",   &fields::Grids::by)
    .def_readwrite("bz",   &fields::Grids::bz)
    .def_readwrite("jx",   &fields::Grids::jx)
    .def_readwrite("jy",   &fields::Grids::jy)
    .def_readwrite("jz",   &fields::Grids::jz)
    .def_readwrite("jx1",  &fields::Grids::jx1)
    .def_readwrite("rho",  &fields::Grids::rho);


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


  py::class_<vlv::VlasovTile, 
             fields::PlasmaTile,
             corgi::Tile, 
             std::shared_ptr<vlv::VlasovTile>
             >(m, "VlasovTile")
    .def(py::init<size_t, size_t, int, size_t, size_t, size_t, size_t>())
    .def_readwrite("dt",     &vlv::VlasovTile::dt)
    .def_readwrite("dx",     &vlv::VlasovTile::dx)
    .def("get_plasma_species", [](vlv::VlasovTile& tile, size_t i, size_t s) 
        { return tile.steps.get(i).at(s); }, py::return_value_policy::reference)
    .def("insert_initial_species", [](vlv::VlasovTile& c, 
                                  std::vector<vlv::PlasmaBlock> species){
        // push twice to initialize both time steps (current and future)
        c.steps.push_back(species);
        c.steps.push_back(species);

        c.Nspecies = species.size();

        })

    .def("clip",         &vlv::VlasovTile::clip)
    .def("cycle",        &vlv::VlasovTile::cycle);



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
    .def("__setitem__", [](AM3d &s, py::tuple indx, float_m v) 
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
      
  //m.def("deriv", &toolbox::deriv<float_m, 3>);
  //m.def("grad",  &toolbox::grad<float_m, 3>);
  //m.def("interp_linear", &toolbox::interp_linear<float_m,3>);
  //m.def("interp_cubic", &toolbox::interp_cubic<float_m,3>);


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
  py::class_<vlv::MomentumSolver<float_m,3>, PyMomentumSolver > vvsol(m, "MomentumSolver");
  vvsol
    .def(py::init<>())
    .def("solve",     &vlv::MomentumSolver<float_m,3>::solve)
    .def("solve_mesh", &vlv::MomentumSolver<float_m,3>::solve_mesh);

  // AMR Lagrangian solver
  py::class_<vlv::AmrMomentumLagrangianSolver<float_m,3>>(m, "AmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());

  py::class_<vlv::GravityAmrMomentumLagrangianSolver<float_m,3>>(m, "GravityAmrMomentumLagrangianSolver", vvsol)
     .def(py::init<>());



  // general interface for spatial solvers
  py::class_<vlv::SpatialSolver<float_m>, PySpatialSolver> vssol(m, "SpatialSolver");
  vssol
    .def(py::init<>())
    .def("solve", &vlv::SpatialSolver<float_m>::solve);


  // AMR Lagrangian solver
  py::class_<vlv::AmrSpatialLagrangianSolver<float_m>>(m, "AmrSpatialLagrangianSolver", vssol)
    .def(py::init<>());


  /// Vlasov tile analyzator
  py::class_<vlv::Analyzator<float_m> >(m, "Analyzator")
    .def(py::init<>());


  declare_mesh<float_m, 0>(m, std::string("Mesh0") );
  declare_mesh<float_m, 1>(m, std::string("Mesh1") );
  declare_mesh<float_m, 3>(m, std::string("Mesh3") );



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


    m.def("initial_step_1d", &vlv::initial_step<1>);
    m.def("initial_step", &vlv::initial_step<3>);

    m.def("step_location", &vlv::step_location);


    m.def("step_velocity_1d", &vlv::step_velocity<1>);
    m.def("step_velocity",   &vlv::step_velocity<3>);
    m.def("step_velocity_with_gravity_1d",   &vlv::step_velocity_with_gravity<1>);

    m.def("analyze",      &vlv::analyze);

    m.def("write_grids", &vlv::write_grids);

    m.def("write_analysis", &vlv::write_analysis);

    m.def("write_mesh", &vlv::write_mesh);

}



