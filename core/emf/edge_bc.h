#pragma once

#include <cstdint>
#include <type_traits>

namespace emf {

/// Boundary condition that sets EMF field values in an edge region.
///
/// The edge region is defined by a direction (0=x, 1=y, 2=z), a side
/// (0=left of position, 1=right of position), and a global coordinate
/// position. Each tile computes its own coverage from the global position.
///
/// Field selection (E, B, or J) is controlled at call time via comm_mode.
/// Per-component masks allow setting only specific components (e.g., only
/// Ey/Ez for a conducting boundary).
///
/// Multiple edge BCs are applied in registration order; later overwrites
/// earlier in overlap regions.
///
/// This struct is trivially copyable (POD) for GPU/SIMD kernel capture.
struct edge_bc {
  using value_type = float;

  std::uint8_t direction = 0;  ///< Edge normal: 0=x, 1=y, 2=z
  std::uint8_t side = 0;       ///< 0=left of position, 1=right of position
  value_type position = 0;     ///< Edge location in global coordinates

  value_type Ex = 0, Ey = 0, Ez = 0;
  value_type Bx = 0, By = 0, Bz = 0;
  value_type Jx = 0, Jy = 0, Jz = 0;

  /// Per-component masks (bit 0=x, 1=y, 2=z). Only masked components are set.
  std::uint8_t E_components = 0b111;
  std::uint8_t B_components = 0b111;
  std::uint8_t J_components = 0b111;
};

static_assert(std::is_trivially_copyable_v<edge_bc>);

}  // namespace emf
