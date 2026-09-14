// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/emf.h"

#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "runko/actions/args.h"
#include "runko/comm/cartesian_grid.h"
#include "runko/coords.h"
#include "runko/emf/common.h"
#include "runko/emf/stencil_coefficients.h"
#include "runko/emf/yee_lattice.h"
#include "runko/simulation_context.h"
#include "runko/tools/config_parser.h"
#include "runko/tools/math.h"
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


void
  register_antenna(runko::simulation_context& sim, emf::antenna_mode mode)
{
  if(not sim.tiles.ctx().contains<emf::antennas>()) {
    sim.tiles.ctx().emplace<emf::antennas>();
  }

  // As the data storage for lap_coeffs in antennas are in std::vector and
  // we use them in order, we can make this little bit nicer by reversing the data
  // and using items from the back:
  if(mode.lap_coeffs) { std::ranges::reverse(mode.lap_coeffs.value()); }

  sim.tiles.ctx().get<emf::antennas>().modes.push_back(mode);
}

/// If modes contain lap coeffs, uses the latest one.
void
  deposit_antenna_current(
    emf::YeeLattice& yee,
    const std::vector<emf::antenna_mode>& modes,
    const runko::global_coordinates_closure<3>& coords,
    const auto global_coords_mins,
    const auto global_coords_maxs,
    const auto cfl,
    emf::antenna_buffers& buffs)
{
  using vec_list = runko::VecList<emf::YeeLattice::value_type>;

  // Fake complex numbers with arrays.
  using complex_list = runko::ScalarList<std::array<emf::YeeLattice::value_type, 2>>;

  const auto num_of_modes = modes.size();
  auto A                  = vec_list(num_of_modes);
  auto K                  = vec_list(num_of_modes);
  auto lap_coeffs         = complex_list(num_of_modes);

  const auto sA_mds          = A.staging_mds();
  const auto sk_mds          = K.staging_mds();
  const auto slap_coeffs_mds = lap_coeffs.staging_mds();

  auto get_wave_vector = [&](const emf::antenna_mode& wm) {
    auto handle_wave_data =
      [&](auto&& data) -> toolbox::Vec3<emf::antenna_mode::value_type> {
      using T = std::decay_t<decltype(data)>;
      if constexpr(std::is_same_v<T, emf::antenna_mode::wave_vector>) {
        return data.k;
      } else if constexpr(std::is_same_v<T, emf::antenna_mode::wave_number>) {
        using F         = emf::antenna_mode::value_type;
        const auto mins = toolbox::Vec3<F>(global_coords_mins);
        const auto maxs = toolbox::Vec3<F>(global_coords_maxs);
        const auto L    = maxs - mins;

        // toolbox::VecD does not have element wise divide.
        const auto tmp = 2 * std::numbers::pi_v<F> * data.n;
        return toolbox::Vec3<F>(tmp[0] / L[0], tmp[1] / L[1], tmp[2] / L[2]);
      }
    };

    return std::visit(handle_wave_data, wm.wave_data);
  };

  for(const auto n: std::views::iota(0uz, num_of_modes)) {
    const auto wave_vector = get_wave_vector(modes[n]);
    for(const auto i: std::views::iota(0uz, 3uz)) {
      sA_mds[n][i] = static_cast<emf::YeeLattice::value_type>(modes[n].A[i]);
      sk_mds[n][i] = static_cast<emf::YeeLattice::value_type>(wave_vector[i]);
    }

    if(modes[n].lap_coeffs and modes[n].lap_coeffs.value().empty()) {
      throw std::logic_error {
        "Can not deposit antenna current, antenna_mode ran out of lap_coeffs!"
      };
    } else if(modes[n].lap_coeffs) {
      const auto z = modes[n].lap_coeffs.value().back();

      slap_coeffs_mds[n][][0] = static_cast<emf::YeeLattice::value_type>(z.real());
      slap_coeffs_mds[n][][1] = static_cast<emf::YeeLattice::value_type>(z.imag());
    } else {
      slap_coeffs_mds[n][][0] = 1;
      slap_coeffs_mds[n][][1] = 0;
    }
  }

  const tyvi::mdgrid_work w {};
  w.sync_from_staging(A).sync_from_staging(K).sync_from_staging(lap_coeffs);

  auto& vec_pot = buffs.vec_pot;
  vec_pot.invalidating_resize(yee.extents_with_halo());

  const auto A_mds          = A.mds();
  const auto K_mds          = K.mds();
  const auto lap_coeffs_mds = lap_coeffs.mds();
  const auto vec_pot_mds    = vec_pot.mds();

  w.for_each_index(vec_pot_mds, [=](const auto idx) {
    // For each is over all indices (including halo region)
    // but the global coordinate map is defined s.t. (0, 0, 0) is located at corner of
    // non-halo region.
    const auto i = static_cast<double>(idx[0]) - emf::halo_size;
    const auto j = static_cast<double>(idx[1]) - emf::halo_size;
    const auto k = static_cast<double>(idx[2]) - emf::halo_size;

    using vec        = toolbox::Vec3<double>;
    const auto x_loc = vec(coords(i + 0.5, j, k));
    const auto y_loc = vec(coords(i, j + 0.5, k));
    const auto z_loc = vec(coords(i, j, k + 0.5));

    vec_pot_mds[idx][0] = 0;
    vec_pot_mds[idx][1] = 0;
    vec_pot_mds[idx][2] = 0;

    for(auto n = 0uz; n < num_of_modes; ++n) {
      // std::complex is contexpr only in c++26 and thus not usable in kernel.
      // Here we do manual complex arithmeitc as a work around.

      const auto phi_x =
        static_cast<emf::YeeLattice::value_type>(toolbox::dot(x_loc, vec(K_mds[n])));
      const auto phi_y =
        static_cast<emf::YeeLattice::value_type>(toolbox::dot(y_loc, vec(K_mds[n])));
      const auto phi_z =
        static_cast<emf::YeeLattice::value_type>(toolbox::dot(z_loc, vec(K_mds[n])));

      const auto x_re = sstd::cos(phi_x);
      const auto x_im = sstd::sin(phi_x);
      const auto y_re = sstd::cos(phi_y);
      const auto y_im = sstd::sin(phi_y);
      const auto z_re = sstd::cos(phi_z);
      const auto z_im = sstd::sin(phi_z);

      const auto w    = std::array<emf::YeeLattice::value_type, 2>(lap_coeffs_mds[n][]);
      const auto w_re = w[0];
      const auto w_im = w[1];


      vec_pot_mds[idx][0] =
        vec_pot_mds[idx][0] + A_mds[n][0] * (w_re * x_re - w_im * x_im);
      vec_pot_mds[idx][1] =
        vec_pot_mds[idx][1] + A_mds[n][1] * (w_re * y_re - w_im * y_im);
      vec_pot_mds[idx][2] =
        vec_pot_mds[idx][2] + A_mds[n][2] * (w_re * z_re - w_im * z_im);
    }
  });

  auto& generated_B = buffs.generated_B;
  generated_B.invalidating_resize(yee.extents_with_halo());

  const auto B_mds = generated_B.mds();

  // We have to calculate B = curl(vec_pot) only in non-halo region + one deep shell in
  // halo region.

  const auto h = emf::halo_size;

  const auto [ex, ey, ez] = yee.extents_wout_halo();
  const auto i1           = std::tuple { h - 1uz, h + ex + 1uz };
  const auto j1           = std::tuple { h - 1uz, h + ey + 1uz };
  const auto k1           = std::tuple { h - 1uz, h + ez + 1uz };

  const auto i1p1 = std::tuple { h, h + ex + 2uz };
  const auto j1p1 = std::tuple { h, h + ey + 2uz };
  const auto k1p1 = std::tuple { h, h + ez + 2uz };

  // FIXME: unify this with emf fdtd2.
  auto curl = [&](
                const auto coeff_,
                const auto out,
                const auto X,
                const auto Xip1,
                const auto Xjp1,
                const auto Xkp1) {
    const auto coeff = static_cast<decltype(out)::element_type::element_type>(coeff_);
    w.for_each_index(out, [=](const auto idx) {
      const auto Dk = Xkp1[idx][1] - X[idx][1];
      const auto Dj = Xjp1[idx][2] - X[idx][2];
      out[idx][0]   = coeff * (Dj - Dk);
    });
    w.for_each_index(out, [=](const auto idx) {
      const auto Di = Xip1[idx][2] - X[idx][2];
      const auto Dk = Xkp1[idx][0] - X[idx][0];
      out[idx][1]   = coeff * (Dk - Di);
    });
    w.for_each_index(out, [=](const auto idx) {
      const auto Dj = Xjp1[idx][0] - X[idx][0];
      const auto Di = Xip1[idx][1] - X[idx][1];
      out[idx][2]   = coeff * (Di - Dj);
    });
  };

  curl(
    1,
    std::submdspan(B_mds, i1, j1, k1),
    std::submdspan(vec_pot_mds, i1, j1, k1),
    std::submdspan(vec_pot_mds, i1p1, j1, k1),
    std::submdspan(vec_pot_mds, i1, j1p1, k1),
    std::submdspan(vec_pot_mds, i1, j1, k1p1));


  // Now curl(B) in non-halo region.
  // We can reuse vec_pot container.

  const auto i   = std::tuple { h, h + ex };
  const auto j   = std::tuple { h, h + ey };
  const auto k   = std::tuple { h, h + ez };
  const auto im1 = std::tuple { h - 1uz, h + ex - 1uz };
  const auto jm1 = std::tuple { h - 1uz, h + ey - 1uz };
  const auto km1 = std::tuple { h - 1uz, h + ez - 1uz };

  // Here we negate cfl, as we put (im1, jm1, km1) instead of (ip1, jp1, kp1).
  curl(
    -cfl,
    std::submdspan(vec_pot_mds, i, j, k),
    std::submdspan(B_mds, i, j, k),
    std::submdspan(B_mds, im1, j, k),
    std::submdspan(B_mds, i, jm1, k),
    std::submdspan(B_mds, i, j, km1));

  yee.deposit_current(w, vec_pot);
  w.wait();
}

ta::sexpr_sender
  deposit_antenna_current(runko::simulation_context& x)
{
  if(not x.tiles.ctx().contains<emf::antenna_buffers>()) {
    x.tiles.ctx().emplace<emf::antenna_buffers>();
  }
  return te::just(std::ref(x)) | te::then([](runko::simulation_context& sim) {
           if(not sim.tiles.ctx().contains<emf::antennas>()) { return ta::null; }
           auto& antennas          = sim.tiles.ctx().get<emf::antennas>();
           auto& buffs             = sim.tiles.ctx().get<emf::antenna_buffers>();
           const auto cfl          = sim.config.template get_or_throw<double>("cfl");
           const auto [mins, maxs] = runko::global_coordinate_extents<3>(sim);
           for(auto&& [_, yee, idx]: sim.view_tiles<
                                     emf::YeeLattice,
                                     runko::cartesian_index<3>,
                                     runko::local_tile_tag>()) {
             const auto coords =
               runko::global_coordinates(sim, idx.template as<double>().data);
             deposit_antenna_current(
               yee,
               antennas.modes,
               coords,
               mins,
               maxs,
               cfl,
               buffs);
           }


           auto pop_lap_coeff = [](emf::antenna_mode& mode) {
             if(mode.lap_coeffs) { mode.lap_coeffs.value().pop_back(); }
           };

           std::ranges::for_each(antennas.modes, pop_lap_coeff);

           return ta::null;
         });
}

constexpr emf::FieldPropagator
  parse_field_propagator(const toolbox::ConfigParser& config)
{

  if(const auto x = config.get<std::string>("field_propagator")) {
    if(x.value() == "fdtd2") {
      return emf::FieldPropagator::fdtd2;
    } else if(x.value() == "stencil") {
      return emf::FieldPropagator::stencil;
    } else {
      std::runtime_error { std::format("unregonized field_propagator: {}", x.value()) };
    }
  }
  throw std::runtime_error {
    "error: configuration parameter missing: field_propagator"
  };
}

constexpr emf::CurrentFilter
  parse_current_filter(const toolbox::ConfigParser& config)
{
  if(const auto x = config.get<std::string>("current_filter")) {
    if(x.value() == "binomial2") {
      return emf::CurrentFilter::binomial2;
    } else if(x.value() == "binomial2_unrolled") {
      return emf::CurrentFilter::binomial2_unrolled;
    } else {
      std::runtime_error { std::format("unregonized current_filter: {}", x.value()) };
    }
  }
  throw std::runtime_error { "configuration parameter missing: current_filter" };
}

constexpr emf::StencilCoeffs
  parse_stencil_coeffs(const toolbox::ConfigParser& config)
{
  // Read stencil coefficients from config.
  // For each coefficient, try per-axis key (stencil_x_name) first, then isotropic
  // (stencil_name). Default to 0.

  auto read_coeff = [&](const std::string& axis_prefix, const std::string& name) {
    // Per-axis key: stencil_x_delta, stencil_y_delta, ...
    if(const auto v = config.get<double>(axis_prefix + name)) {
      return static_cast<float>(v.value());
    }
    // Isotropic key: stencil_delta, stencil_gamma, ...
    if(const auto v = config.get<double>("stencil_" + name)) {
      return static_cast<float>(v.value());
    }
    return 0.0f;
  };

  auto read_axis = [&](const std::string& axis_prefix) -> emf::StencilAxisCoeffs {
    emf::StencilAxisCoeffs c {};
    c.M[1][0] = read_coeff(axis_prefix, "delta");
    c.M[2][0] = read_coeff(axis_prefix, "gamma");
    c.M[0][1] = read_coeff(axis_prefix, "beta_p1");
    c.M[0][2] = read_coeff(axis_prefix, "beta_p2");
    c.M[1][1] = read_coeff(axis_prefix, "beta2_p1");
    c.M[1][2] = read_coeff(axis_prefix, "beta2_p2");
    c.M[2][1] = read_coeff(axis_prefix, "beta3_p1");
    c.M[2][2] = read_coeff(axis_prefix, "beta3_p2");
    c.M[0][3] = read_coeff(axis_prefix, "zeta_p1");
    c.M[0][4] = read_coeff(axis_prefix, "zeta_p2");
    c.M[1][3] = read_coeff(axis_prefix, "zeta2_p1");
    c.M[1][4] = read_coeff(axis_prefix, "zeta2_p2");
    c.M[2][3] = read_coeff(axis_prefix, "zeta3_p1");
    c.M[2][4] = read_coeff(axis_prefix, "zeta3_p2");
    // Set alpha from normalization
    c.M[0][0] = c.alpha();
    return c;
  };

  emf::StencilCoeffs coeffs {};
  coeffs.axis[0] = read_axis("stencil_x_");
  coeffs.axis[1] = read_axis("stencil_y_");
  coeffs.axis[2] = read_axis("stencil_z_");

  return coeffs;
}


ta::sexpr_sender
  push_e(runko::simulation_context& x)
{

  const auto cfl  = x.config.template get_or_throw<double>("cfl");
  const auto prop = x.get_n_set_config<emf::FieldPropagator>(&parse_field_propagator);

  using vt = emf::YeeLattice::value_type;

  return te::just(std::ref(x)) | te::then([prop, cfl](runko::simulation_context& sim) {
           for(auto&& [_, yee]:
               sim.view_tiles<emf::YeeLattice, runko::local_tile_tag>()) {
             switch(prop) {
               case FieldPropagator::fdtd2:
               case FieldPropagator::stencil:
                 yee.push_e_fdtd2(static_cast<vt>(cfl));
                 break;
               default:
                 throw std::logic_error {
                   "internal error: unregonized FieldPropagator"
                 };
             }
           }

           return ta::null;
         });
}


ta::sexpr_sender
  push_half_b(runko::simulation_context& x)
{

  const auto cfl  = x.config.template get_or_throw<double>("cfl");
  const auto prop = x.get_n_set_config<emf::FieldPropagator>(&parse_field_propagator);

  if(prop == emf::FieldPropagator::stencil) {
    std::ignore = x.set_config<emf::StencilCoeffs>(&parse_stencil_coeffs);
  }

  using vt = emf::YeeLattice::value_type;

  return te::just(std::ref(x)) | te::then([prop, cfl](runko::simulation_context& sim) {
           for(auto&& [_, yee]:
               sim.view_tiles<emf::YeeLattice, runko::local_tile_tag>()) {
             switch(prop) {
               case FieldPropagator::fdtd2:
                 yee.push_b_fdtd2(static_cast<vt>(cfl / 2));
                 break;
               case FieldPropagator::stencil:
                 yee.push_b_stencil(
                   static_cast<vt>(cfl / 2),
                   sim.get_config<emf::StencilCoeffs>());
                 break;
               default:
                 throw std::logic_error {
                   "internal error: unregonized FieldPropagator"
                 };
             }
           }

           return ta::null;
         });
}

ta::sexpr_sender
  filter_current(runko::simulation_context& x)
{
  const auto filter = x.get_n_set_config<emf::CurrentFilter>(&parse_current_filter);

  return te::just(std::ref(x)) | te::then([filter](runko::simulation_context& sim) {
           for(auto&& [_, yee]:
               sim.view_tiles<emf::YeeLattice, runko::local_tile_tag>()) {
             switch(filter) {
               case emf::CurrentFilter::binomial2:
                 yee.filter_current_binomial2();
                 break;
               case emf::CurrentFilter::binomial2_unrolled:
                 yee.filter_current_binomial2_unrolled();
                 break;
               default:
                 throw std::logic_error {
                   "filter_current internal error: unregonized current filter."
                 };
             }
           }

           return ta::null;
         });
}

}  // namespace emf
