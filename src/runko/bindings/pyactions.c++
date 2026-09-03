#include "corgi/corgi.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/actions/env.h"
#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/runtime.h"
#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"
#include "tyvi/actions_eval.h"

#include <chrono>
#include <exception>
#include <pika/execution.hpp>
#include <pika/init.hpp>
#include <pika/thread.hpp>
#include <print>
#include <ranges>
#include <thread>
#include <tuple>
#include <vector>

namespace {

namespace py = pybind11;
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace rn = std::ranges;
namespace rv = std::views;

ta::sexpr
  parse_element(const py::handle &obj)
{
  if(py::isinstance<py::int_>(obj)) {
    return obj.cast<long>();
  } else if(py::isinstance<py::str>(obj)) {
    return obj.cast<std::string>();
  } else if(py::isinstance<ta::intrinsic>(obj)) {
    return obj.cast<ta::intrinsic>();
  } else if(py::isinstance<runko::symbol>(obj)) {
    return obj.cast<runko::symbol>();
  } else if(py::isinstance<runko::comm_mode>(obj)) {
    return obj.cast<runko::comm_mode>();
  } else if(py::isinstance<py::tuple>(obj)) {
    const auto tup = obj.cast<py::tuple>();

    if(rn::empty(tup)) { return ta::cons(); }

    const auto n = rn::size(tup);
    auto tail    = ta::sexpr { ta::cons(parse_element(tup[n - 1uz]), ta::null) };
    for(const auto i: rv::iota(0uz, n) | rv::reverse | rv::drop(1)) {
      tail = ta::cons(parse_element(tup[i]), std::move(tail));
    }

    return tail;
  } else {
    throw std::runtime_error { "Trying to parse unsupported type." };
  }
}

void
  empty_context_eval(const py::handle &body_py)
{
  const auto body = parse_element(body_py);

  [[maybe_unused]] runko::RuntimeActivator _;

  try {
    tyvi::this_thread::sync_wait(ta::eval<runko::symbol>(body, runko::build_stdenv()));
  } catch(const std::exception &e) {
    std::println("Evaluation exception in empty_context_eval: {}", e.what());
  }
}

void
  simulation_context_eval(const py::handle &body_py, runko::simulation_context &sim)
{
  const auto body = parse_element(body_py);

  [[maybe_unused]] runko::RuntimeActivator _;

  try {
    tyvi::this_thread::sync_wait(ta::eval<runko::symbol>(body, runko::sim_env(sim)));
  } catch(const std::exception &e) {
    std::println("Evaluation exception in simulation_context_eval: {}", e.what());
  }
}

void
  add_init_tiles(runko::simulation_context &sim, const py::handle &tiles)
{
  for(const auto &tile: tiles) {
    const auto idx =
      tile.attr("_idx").cast<std::tuple<std::size_t, std::size_t, std::size_t>>();

    const auto id    = sim.tiles.create();
    using index_type = runko::cartesian_index<3>;
    sim.tiles
      .emplace<index_type>(id, std::get<0>(idx), std::get<1>(idx), std::get<2>(idx));
    sim.tiles.emplace<runko::local_tile_tag>(id);

    const auto extents =
      tile.attr("_extents_wout_halo").cast<std::vector<std::ptrdiff_t>>();
    auto yee = emf::YeeLattice(emf::make_yee_ctor_args_from_vector(extents));

    const auto Ex = tile.attr("_Ex").cast<py::array_t<double>>();
    const auto Ey = tile.attr("_Ey").cast<py::array_t<double>>();
    const auto Ez = tile.attr("_Ez").cast<py::array_t<double>>();
    const auto Bx = tile.attr("_Bx").cast<py::array_t<double>>();
    const auto By = tile.attr("_By").cast<py::array_t<double>>();
    const auto Bz = tile.attr("_Bz").cast<py::array_t<double>>();
    const auto Jx = tile.attr("_Jx").cast<py::array_t<double>>();
    const auto Jy = tile.attr("_Jy").cast<py::array_t<double>>();
    const auto Jz = tile.attr("_Jz").cast<py::array_t<double>>();

    const auto Exv = Ex.template unchecked<3>();
    const auto Eyv = Ey.template unchecked<3>();
    const auto Ezv = Ez.template unchecked<3>();
    const auto Bxv = Bx.template unchecked<3>();
    const auto Byv = By.template unchecked<3>();
    const auto Bzv = Bz.template unchecked<3>();
    const auto Jxv = Jx.template unchecked<3>();
    const auto Jyv = Jy.template unchecked<3>();
    const auto Jzv = Jz.template unchecked<3>();

    yee.set_EBJ([&](const std::size_t i, const std::size_t j, const std::size_t k) {
      return emf::YeeLatticeFieldsAtPoint { .Ex = Exv(i, j, k),
                                            .Ey = Eyv(i, j, k),
                                            .Ez = Ezv(i, j, k),
                                            .Bx = Bxv(i, j, k),
                                            .By = Byv(i, j, k),
                                            .Bz = Bzv(i, j, k),
                                            .Jx = Jxv(i, j, k),
                                            .Jy = Jyv(i, j, k),
                                            .Jz = Jzv(i, j, k) };
    });


    sim.tiles.emplace<emf::YeeLattice>(id, std::move(yee));
  }

  runko::set_cartesian_neighbors<3>(sim);
  runko::set_cartesian_comm_infos<3>(sim);
}


/// MVP implementation.
auto
  to_ndarrays(const emf::YeeLattice::YeeLatticeHostCopy &lattice)
{

  const auto grid_shape = std::array { lattice.grid_extents().extent(0),
                                       lattice.grid_extents().extent(1),
                                       lattice.grid_extents().extent(2) };

  auto Ex = py::array_t<double, py::array::c_style>(grid_shape);
  auto Ey = py::array_t<double, py::array::c_style>(grid_shape);
  auto Ez = py::array_t<double, py::array::c_style>(grid_shape);
  auto Bx = py::array_t<double, py::array::c_style>(grid_shape);
  auto By = py::array_t<double, py::array::c_style>(grid_shape);
  auto Bz = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jx = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jy = py::array_t<double, py::array::c_style>(grid_shape);
  auto Jz = py::array_t<double, py::array::c_style>(grid_shape);

  auto Exv = Ex.template mutable_unchecked<3>();
  auto Eyv = Ey.template mutable_unchecked<3>();
  auto Ezv = Ez.template mutable_unchecked<3>();
  auto Bxv = Bx.template mutable_unchecked<3>();
  auto Byv = By.template mutable_unchecked<3>();
  auto Bzv = Bz.template mutable_unchecked<3>();
  auto Jxv = Jx.template mutable_unchecked<3>();
  auto Jyv = Jy.template mutable_unchecked<3>();
  auto Jzv = Jz.template mutable_unchecked<3>();

  for(const auto mds = lattice.mds(); const auto idx: tyvi::sstd::index_space(mds)) {
    const auto [i, j, k] = idx;
    const auto F         = mds[idx][];

    Exv(i, j, k) = F.Ex;
    Eyv(i, j, k) = F.Ey;
    Ezv(i, j, k) = F.Ez;
    Bxv(i, j, k) = F.Bx;
    Byv(i, j, k) = F.By;
    Bzv(i, j, k) = F.Bz;
    Jxv(i, j, k) = F.Jx;
    Jyv(i, j, k) = F.Jy;
    Jzv(i, j, k) = F.Jz;
  }

  return std::tuple { std::tuple { std::move(Ex), std::move(Ey), std::move(Ez) },
                      std::tuple { std::move(Bx), std::move(By), std::move(Bz) },
                      std::tuple { std::move(Jx), std::move(Jy), std::move(Jz) } };
}

auto
  get_EBJ_with_halo(
    runko::simulation_context &sim,
    const runko::simulation_context::tile_id_type id)
{
  if(const auto p = sim.tiles.try_get<emf::YeeLattice>(id)) {
    return to_ndarrays(p->get_EBJ_with_halo());
  } else {
    throw std::runtime_error(
      "Trying to invoke get_EBJ_with_halo of a tile which do not have YeeLattice.");
  }
}

auto
  debug_cartesian_grid(const runko::simulation_context &sim)
{
  for(const auto [id, idx, neighs]: sim.view_tiles<
                                    runko::cartesian_index<3>,
                                    runko::cartesian_neighbors<3>,
                                    runko::local_tile_tag>()) {
    std::println("local tile at: {} {} {}", idx[0], idx[1], idx[2]);
    std::println("has neighbors:");

    using neigh_type = runko::grid_neighbor<3>;

    for(const auto dir:
        rv::iota(0uz, std::pow(3uz, 3uz)) |
          rv::transform([&](const auto n) { return neigh_type::from_index(n); }) |
          rv::filter(
            [](const auto &x) { return x != runko::grid_neighbor_origo<3>; })) {

      const auto dir_vec = dir.to_vec<int>();
      std::print("{} {} {}: ", dir_vec[0], dir_vec[1], dir_vec[2]);
      if(const auto neigh_id = neighs.get(dir)) {
        const auto i = sim.tiles.try_get<runko::cartesian_index<3>>(neigh_id.value());
        if(i) {
          std::println("{} {} {}", (*i)[0], (*i)[1], (*i)[2]);
        } else {
          std::println("id from cartesian_neighbors did not point to valid tile");
        }

      } else {
        std::println("no cartesian index");
      }
    }
  }
}

}  // namespace


namespace actions {

void
  bind_actions(py::module &m_sub)
{
  py::enum_<ta::intrinsic>(m_sub, "intrinsic")
    .value("car", ta::intrinsic::car)
    .value("cdr", ta::intrinsic::cdr)
    .value("quote", ta::intrinsic::quote)
    .export_values();

  py::enum_<runko::symbol>(m_sub, "symbol")
    .value("print", runko::symbol::print)
    .value("println", runko::symbol::println)
    .value("version", runko::symbol::version)
    .value("mt_showcase", runko::symbol::mt_showcase)
    .value("comm_local", runko::symbol::comm_local)
    .value("comm_external", runko::symbol::comm_external)
    .value("current_context", runko::symbol::current_context)
    .export_values();

  m_sub.def("empty_context_eval", &::empty_context_eval);

  py::class_<runko::RuntimeInstance>(m_sub, "RuntimeInstance").def(py::init<>());

  py::class_<runko::simulation_context::tile_id_type>(m_sub, "TileID")
    .def(py::init<>());

  py::class_<runko::simulation_context>(m_sub, "SimulationContext")
    .def(py::init([](const py::handle &config) {
      return runko::simulation_context { .config { toolbox::ConfigParser(config) },
                                         .tiles {} };
    }))
    .def(
      "eval",
      [](runko::simulation_context &sim, const py::handle &body_py) {
        simulation_context_eval(body_py, sim);
      })
    .def("add_init_tiles", &add_init_tiles)
    .def(
      "get_local_tile_ids",
      [](const runko::simulation_context &sim) {
        auto ids = std::vector<runko::simulation_context::tile_id_type> {};
        for(const auto [id]: sim.view_tiles<runko::local_tile_tag>()) {
          ids.push_back(id);
        }
        return ids;
      })
    .def("get_EBJ_with_halo", &get_EBJ_with_halo)
    .def("debug_cartesian_grid", &debug_cartesian_grid);
}
}  // namespace actions
