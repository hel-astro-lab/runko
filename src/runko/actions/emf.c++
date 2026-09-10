// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/emf.h"

#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "runko/actions/args.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/coords.h"
#include "runko/emf/common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/config_parser.h"
#include "tyvi/execution.h"
#include "tyvi/mdspan.h"

#include <algorithm>
#include <functional>
#include <print>
#include <tuple>

namespace emf {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace py = pybind11;

void
  ensure_constructed_yee_lattices(
    const std::reference_wrapper<runko::simulation_context> sim)
{
  auto where_to_construct = std::vector<runko::simulation_context::tile_id_type> {};

  for(auto&& id:
      sim.get().tiles.view<runko::local_tile_tag>(entt::exclude<emf::YeeLattice>)) {
    where_to_construct.push_back(id);
  }
  if(where_to_construct.empty()) { return; }

  const auto n_cells =
    toolbox::get_extent_list(sim.get().config, "n_cells_per_tile", 3);
  const auto args = emf::make_yee_ctor_args_from_vector(n_cells);
  std::ranges::for_each(where_to_construct, [&](const auto x) {
    sim.get().tiles.template emplace<emf::YeeLattice>(x, args);
  });
}

ta::sexpr_sender
  set_EBJ(
    const std::reference_wrapper<runko::simulation_context> sim,
    vector_field_function E,
    vector_field_function B,
    vector_field_function J)
{
  return te::just() | te::then([=] {
           for(auto&& [id, yee, idx]: sim.get()
                                        .template view_tiles<
                                          emf::YeeLattice,
                                          const runko::cartesian_index<3>,
                                          runko::local_tile_tag>()) {
             const auto gc =
               runko::global_coordinates(sim.get(), idx.template as<double>().data);
             auto f =
               [&](const std::size_t i, const std::size_t j, const std::size_t k) {
                 const auto [x, y, z] = gc(i, j, k);
                 const auto ex        = E(x + 0.5, y, z);
                 const auto ey        = E(x, y + 0.5, z);
                 const auto ez        = E(x, y, z + 0.5);
                 const auto bx        = B(x, y + 0.5, z + 0.5);
                 const auto by        = B(x + 0.5, y, z + 0.5);
                 const auto bz        = B(x + 0.5, y + 0.5, z);
                 const auto jx        = J(x + 0.5, y, z);
                 const auto jy        = J(x, y + 0.5, z);
                 const auto jz        = J(x, y, z + 0.5);

                 return emf::YeeLatticeFieldsAtPoint { .Ex = std::get<0>(ex),
                                                       .Ey = std::get<1>(ey),
                                                       .Ez = std::get<2>(ez),
                                                       .Bx = std::get<0>(bx),
                                                       .By = std::get<1>(by),
                                                       .Bz = std::get<2>(bz),
                                                       .Jx = std::get<0>(jx),
                                                       .Jy = std::get<1>(jy),
                                                       .Jz = std::get<2>(jz) };
               };
             yee.set_EBJ(f);
           }
           return ta::null;
         });
}

ta::sexpr_sender
  batch_set_EBJ(
    const std::reference_wrapper<runko::simulation_context> sim,
    batch_vector_field_function Ex,
    batch_vector_field_function Ey,
    batch_vector_field_function Ez,
    batch_vector_field_function Bx,
    batch_vector_field_function By,
    batch_vector_field_function Bz,
    batch_vector_field_function Jx,
    batch_vector_field_function Jy,
    batch_vector_field_function Jz)
{
  return te::just() | te::then([=] {
           for(auto&& [id, yee, idx]: sim.get()
                                        .template view_tiles<
                                          emf::YeeLattice,
                                          const runko::cartesian_index<3>,
                                          runko::local_tile_tag>()) {
             const auto global_coordinates =
               runko::global_coordinates(sim.get(), idx.template as<double>().data);

             const auto e = yee.extents_wout_halo();

             auto x   = pybind11::array_t<double>(e);
             auto y   = pybind11::array_t<double>(e);
             auto z   = pybind11::array_t<double>(e);
             auto xp5 = pybind11::array_t<double>(e);
             auto yp5 = pybind11::array_t<double>(e);
             auto zp5 = pybind11::array_t<double>(e);

             auto xv   = x.template mutable_unchecked<3>();
             auto yv   = y.template mutable_unchecked<3>();
             auto zv   = z.template mutable_unchecked<3>();
             auto xp5v = xp5.template mutable_unchecked<3>();
             auto yp5v = yp5.template mutable_unchecked<3>();
             auto zp5v = zp5.template mutable_unchecked<3>();


             for(const auto [i, j, k]: tyvi::sstd::index_space(
                   std::mdspan((int*)nullptr, e[0], e[1], e[2]))) {

               const auto ii           = static_cast<double>(i);
               const auto jj           = static_cast<double>(j);
               const auto kk           = static_cast<double>(k);
               const auto [gx, gy, gz] = global_coordinates(i, j, k);
               const auto [gxp5, gyp5, gzp5] =
                 global_coordinates(ii + 0.5, jj + 0.5, kk + 0.5);
               xv(i, j, k)   = gx;
               yv(i, j, k)   = gy;
               zv(i, j, k)   = gz;
               xp5v(i, j, k) = gxp5;
               yp5v(i, j, k) = gyp5;
               zp5v(i, j, k) = gzp5;
             }

             const auto ex = Ex(xp5, y, z);
             const auto ey = Ey(x, yp5, z);
             const auto ez = Ez(x, y, zp5);

             const auto bx = Bx(x, yp5, zp5);
             const auto by = By(xp5, y, zp5);
             const auto bz = Bz(xp5, yp5, z);

             const auto jx = Jx(xp5, y, z);
             const auto jy = Jy(x, yp5, z);
             const auto jz = Jz(x, y, zp5);

             const auto assert_shape = [&](const auto& A) {
               if(
                 static_cast<std::size_t>(A.shape(0)) != e[0] or
                 static_cast<std::size_t>(A.shape(1)) != e[1] or
                 static_cast<std::size_t>(A.shape(2)) != e[2]) {
                 throw std::runtime_error {
                   "Batch field setter returned array with incorrect shape!"
                 };
               }
             };

             assert_shape(ex);
             assert_shape(ey);
             assert_shape(ez);
             assert_shape(bx);
             assert_shape(by);
             assert_shape(bz);
             assert_shape(jx);
             assert_shape(jy);
             assert_shape(jz);

             const auto exv = ex.template unchecked<3>();
             const auto eyv = ey.template unchecked<3>();
             const auto ezv = ez.template unchecked<3>();
             const auto bxv = bx.template unchecked<3>();
             const auto byv = by.template unchecked<3>();
             const auto bzv = bz.template unchecked<3>();
             const auto jxv = jx.template unchecked<3>();
             const auto jyv = jy.template unchecked<3>();
             const auto jzv = jz.template unchecked<3>();


             auto f =
               [&](const std::size_t i, const std::size_t j, const std::size_t k) {
                 return YeeLatticeFieldsAtPoint { .Ex = exv(i, j, k),
                                                  .Ey = eyv(i, j, k),
                                                  .Ez = ezv(i, j, k),
                                                  .Bx = bxv(i, j, k),
                                                  .By = byv(i, j, k),
                                                  .Bz = bzv(i, j, k),
                                                  .Jx = jxv(i, j, k),
                                                  .Jy = jyv(i, j, k),
                                                  .Jz = jzv(i, j, k) };
               };

             yee.set_EBJ(f);
           }
           return ta::null;
         });
}


ta::sexpr_sender
  add_current(const std::reference_wrapper<runko::simulation_context> sim)
{
  return te::just() | te::then([=] {
           for(auto&& [id, yee]:
               sim.get()
                 .template view_tiles<emf::YeeLattice, runko::local_tile_tag>()) {
             yee.add_current();
           }
           return ta::null;
         });
}

}  // namespace emf
