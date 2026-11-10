// Copyright 2016 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/communication_common.h"
#include "runko/emf/antenna.h"
#include "runko/emf/edge_bc.h"
#include "runko/tools/config_parser.h"
#include "tyvi/mdgrid_buffer.h"
#include "tyvi/mdspan.h"

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

//--------------------------------------------------

namespace emf {

namespace py = pybind11;

void
  bind_emf(py::module& m_sub)
{

  //--------------------------------------------------
  // 1D bindings
  // TODO

  //--------------------------------------------------
  // 2D bindings
  // TODO

  //--------------------------------------------------
  // 3D bindings
  py::module m_3d = m_sub.def_submodule("threeD", "3D specializations");

  using antenna_pyvec = pybind11::array_t<emf::antenna_mode::value_type>;
  using antenna_complex_pyvec =
    pybind11::array_t<std::complex<emf::antenna_mode::value_type>>;
  py::class_<emf::antenna_mode>(m_3d, "antenna_mode")
    .def(
      py::init([](
                 antenna_pyvec A,
                 std::optional<antenna_pyvec> k,
                 std::optional<antenna_pyvec> n,
                 std::optional<antenna_complex_pyvec> lap_coeffs) {
        if((k and n) or (not k and not n)) {
          throw std::runtime_error(
            "antenna_mode expects k or n to be defined but not both.");
        }


        auto assert_3d_vec = [](auto& x) {
          if(x.ndim() != 1) {
            throw std::runtime_error(
              "Antenna expects A and k/n to be rank-1 arrays (specifically 3D "
              "vectors).");
          }

          if(x.shape(0) != 3) {
            throw std::runtime_error("Antenna expects A and k/n to be 3D vectors.");
          }
        };

        const auto wave_data =
          std::invoke([&] -> decltype(emf::antenna_mode::wave_data) {
            if(k) {
              assert_3d_vec(k.value());
              const auto kv = k.value().template unchecked<1>();
              return emf::antenna_mode::wave_vector { { kv(0), kv(1), kv(2) } };
            } else {
              assert_3d_vec(n.value());
              const auto nv = n.value().template unchecked<1>();
              return emf::antenna_mode::wave_number { { nv(0), nv(1), nv(2) } };
            }
          });

        assert_3d_vec(A);
        const auto Av = A.template unchecked<1>();

        auto to_stdvec = [](antenna_complex_pyvec& p)
          -> std::optional<std::vector<std::complex<emf::antenna_mode::value_type>>> {
          if(p.ndim() != 1) {
            throw std::runtime_error { "lap_coeffs must be 1D array." };
          }

          const auto N  = static_cast<std::size_t>(p.shape(0));
          auto vec      = std::vector<std::complex<emf::antenna_mode::value_type>>(N);
          const auto pv = p.template unchecked<1>();

          for(auto i = 0uz; i < N; ++i) { vec[i] = pv(i); }
          return vec;
        };

        return emf::antenna_mode { .A { Av(0), Av(1), Av(2) },
                                   .wave_data { wave_data },
                                   .lap_coeffs { lap_coeffs.and_then(to_stdvec) } };
      }),
      py::kw_only(),
      py::arg("A"),
      py::arg("k")          = std::optional<antenna_pyvec> {},
      py::arg("n")          = std::optional<antenna_pyvec> {},
      py::arg("lap_coeffs") = std::optional<antenna_complex_pyvec> {});

  // edge boundary condition struct
  using EBC = emf::edge_bc;
  py::class_<EBC>(m_3d, "edge_bc")
    .def(
      py::init([](
                 std::uint8_t direction,
                 std::uint8_t side,
                 EBC::value_type position,
                 EBC::value_type Ex,
                 EBC::value_type Ey,
                 EBC::value_type Ez,
                 EBC::value_type Bx,
                 EBC::value_type By,
                 EBC::value_type Bz,
                 EBC::value_type Jx,
                 EBC::value_type Jy,
                 EBC::value_type Jz,
                 std::uint8_t E_components,
                 std::uint8_t B_components,
                 std::uint8_t J_components) {
        return EBC { direction, side, position,     Ex,           Ey,
                     Ez,        Bx,   By,           Bz,           Jx,
                     Jy,        Jz,   E_components, B_components, J_components };
      }),
      py::kw_only(),
      py::arg("direction")    = std::uint8_t { 0 },
      py::arg("side")         = std::uint8_t { 0 },
      py::arg("position")     = EBC::value_type { 0 },
      py::arg("Ex")           = EBC::value_type { 0 },
      py::arg("Ey")           = EBC::value_type { 0 },
      py::arg("Ez")           = EBC::value_type { 0 },
      py::arg("Bx")           = EBC::value_type { 0 },
      py::arg("By")           = EBC::value_type { 0 },
      py::arg("Bz")           = EBC::value_type { 0 },
      py::arg("Jx")           = EBC::value_type { 0 },
      py::arg("Jy")           = EBC::value_type { 0 },
      py::arg("Jz")           = EBC::value_type { 0 },
      py::arg("E_components") = std::uint8_t { 0b111 },
      py::arg("B_components") = std::uint8_t { 0b111 },
      py::arg("J_components") = std::uint8_t { 0b111 })
    .def_readwrite("direction", &EBC::direction)
    .def_readwrite("side", &EBC::side)
    .def_readwrite("position", &EBC::position)
    .def_readwrite("Ex", &EBC::Ex)
    .def_readwrite("Ey", &EBC::Ey)
    .def_readwrite("Ez", &EBC::Ez)
    .def_readwrite("Bx", &EBC::Bx)
    .def_readwrite("By", &EBC::By)
    .def_readwrite("Bz", &EBC::Bz)
    .def_readwrite("Jx", &EBC::Jx)
    .def_readwrite("Jy", &EBC::Jy)
    .def_readwrite("Jz", &EBC::Jz)
    .def_readwrite("E_components", &EBC::E_components)
    .def_readwrite("B_components", &EBC::B_components)
    .def_readwrite("J_components", &EBC::J_components);
}
}  // namespace emf
