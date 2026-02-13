#pragma once


#include "tools/vector.h"

#include <complex>
#include <optional>
#include <vector>


namespace emf {

/// Represents single vector potential Fourier mode.
///
/// Values are in code units.
///
/// Resulting current is calculated as:
///
///     J = Re((c / 4pi) * curl(curl(A * phi[n] * exp(i * dot(k, x)))))
///
/// where phi = phase_evolution, n = lap and x = location in global coordinates.
///
/// In code units:
///
///     J = Re(cfl * curl(curl(A * phi[n] * exp(i * dot(k, x)))))
struct antenna_mode {
  using value_type = double;

  toolbox::Vec3<value_type> A {};
  toolbox::Vec3<value_type> k {};
  // std::optional<std::vector<std::complex<value_type>>> phase_evolution{};
};

}  // namespace emf
