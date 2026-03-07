#pragma once

namespace pic {

/// Reflecting wall (conducting piston) in the yz-plane.
///
/// Wall normal is the x-direction.
/// Particles crossing the wall get reflected (u_x -> -u_x for stationary wall).
/// E_y = E_z = 0 behind the wall (conducting boundary condition).
/// Correction currents are deposited to cancel current behind the wall.
struct reflector_wall {
  using value_type = float;

  /// Current x-coordinate of the wall in global (code-unit) coordinates.
  value_type walloc = 0;

  /// Wall velocity: beta = v_wall / c.
  value_type betawall = 0;

  /// Lorentz factor of the wall: gamma = 1/sqrt(1 - beta^2).
  value_type gammawall = 1;
};

}  // namespace pic
