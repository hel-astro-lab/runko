#pragma once

#include <array>

namespace emf {

/// Per-axis coefficient matrix for the extended FDTD stencil (N=2, M=2).
///
/// The 3×5 matrix M maps the Dirichlet kernel basis d(κ) and transverse
/// cosine basis T(κ_perp) to the amplification factor A = d^T · M · T.
///
/// Layout:
///
///          col 0    col 1       col 2       col 3       col 4
///          (axial)  (cos κ_p1)  (cos κ_p2)  (cos 2κ_p1) (cos 2κ_p2)
/// row 0 D₀  alpha    beta_p1     beta_p2     zeta_p1     zeta_p2
/// row 1 D₁  delta    beta2_p1    beta2_p2    zeta2_p1    zeta2_p2
/// row 2 D₂  gamma    beta3_p1    beta3_p2    zeta3_p1    zeta3_p2
///
/// Perpendicular directions follow cyclic ordering:
///   x-axis: perp1=y, perp2=z
///   y-axis: perp1=z, perp2=x
///   z-axis: perp1=x, perp2=y
///
/// M[0][0] (alpha) is auto-computed from the normalization constraint
/// d(0)^T · M · T(0,0) = 1, so the user sets only M[row>0][0] and M[*][col>0].
struct StencilAxisCoeffs {
  float M[3][5] = {};

  /// Compute alpha from the normalization constraint.
  ///
  /// d(0) = (1, 3, 5)^T,  T(0,0) = (1, 2, 2, 2, 2)^T
  /// d(0)^T · M · T(0,0) = 1
  ///
  /// alpha = 1 - 3*delta - 5*gamma
  ///   - 2*(beta_p1 + beta_p2 + zeta_p1 + zeta_p2)
  ///   - 6*(beta2_p1 + beta2_p2 + zeta2_p1 + zeta2_p2)
  ///   - 10*(beta3_p1 + beta3_p2 + zeta3_p1 + zeta3_p2)
  constexpr float alpha() const
  {
    // d(0) weights for rows: (1, 3, 5)
    // T(0,0) weights for cols: (1, 2, 2, 2, 2)
    // Sum over all entries except M[0][0]:
    //   row_n contribution = d(0)[n] * (M[n][0]*1 + sum_{c>0} M[n][c]*2)
    // Solve for M[0][0] from total = 1.

    float sum = 0.0f;

    // Row 1 (delta row): weight 3
    sum += 3.0f * (M[1][0] + 2.0f * (M[1][1] + M[1][2] + M[1][3] + M[1][4]));
    // Row 2 (gamma row): weight 5
    sum += 5.0f * (M[2][0] + 2.0f * (M[2][1] + M[2][2] + M[2][3] + M[2][4]));
    // Row 0 transverse (cols 1-4): weight 1
    sum += 1.0f * 2.0f * (M[0][1] + M[0][2] + M[0][3] + M[0][4]);

    return 1.0f - sum;
  }
};


/// Full stencil coefficients for all three derivative axes.
struct StencilCoeffs {
  std::array<StencilAxisCoeffs, 3> axis {};  // 0=x, 1=y, 2=z

  /// Create isotropic coefficients (all axes identical).
  static constexpr StencilCoeffs isotropic(StencilAxisCoeffs c) { return { { c, c, c } }; }
};

}  // namespace emf
