#include "core/emf/yee_lattice.h"
#include "tyvi/mdgrid.h"

namespace emf {

void
  YeeLattice::push_b_stencil(const value_type dt, const StencilCoeffs& coeffs)
{
  tyvi::mdgrid_work w {};
  push_b_stencil(w, dt, coeffs);
  w.wait();
}

void
  YeeLattice::push_b_stencil(
    const tyvi::mdgrid_work& w,
    const value_type dt,
    const StencilCoeffs& coeffs)
{
  const auto h = static_cast<runko::index_t>(halo_size);

  const auto Bmds      = nonhalo_submds(B_.mds());
  const auto Emds_halo = E_.mds();

  // Per-axis coefficient matrices.
  // For derivative D*_i, axis[i] holds the 3x5 coefficient matrix.
  // Perpendicular directions are cyclic: x->(y,z), y->(z,x), z->(x,y).
  const auto cx = coeffs.axis[0];
  const auto cy = coeffs.axis[1];
  const auto cz = coeffs.axis[2];

  // === Bx: Bx += dt * (D*_z Ey - D*_y Ez) ===
  //
  // D*_z Ey: z is axial direction, perp1=x, perp2=y
  //   axial offsets in k: n=0 -> (k+1, k), n=1 -> (k+2, k-1), n=2 -> (k+3, k-2)
  //   perp1(x) offsets: ±m in i
  //   perp2(y) offsets: ±m in j
  //
  // D*_y Ez: y is axial direction, perp1=z, perp2=x
  //   axial offsets in j: n=0 -> (j+1, j), n=1 -> (j+2, j-1), n=2 -> (j+3, j-2)
  //   perp1(z) offsets: ±m in k
  //   perp2(x) offsets: ±m in i

  w.for_each_index(Bmds, [=](const auto idx) {
    const auto i = idx[0] + h;
    const auto j = idx[1] + h;
    const auto k = idx[2] + h;
    const auto E = Emds_halo;

    // --- D*_z Ey (component 1) using cz coefficients ---
    // z axial, perp1=x (cols 1,3), perp2=y (cols 2,4)
    const auto DzEy =
      // n=0: axial k+1, k
        cz.M[0][0] *  (E[i  , j  , k+1][1] - E[i  , j  , k  ][1])
      + cz.M[0][1] * ((E[i+1, j  , k+1][1] + E[i-1, j  , k+1][1])
                    - (E[i+1, j  , k  ][1] + E[i-1, j  , k  ][1]))
      + cz.M[0][2] * ((E[i  , j+1, k+1][1] + E[i  , j-1, k+1][1])
                    - (E[i  , j+1, k  ][1] + E[i  , j-1, k  ][1]))
      + cz.M[0][3] * ((E[i+2, j  , k+1][1] + E[i-2, j  , k+1][1])
                    - (E[i+2, j  , k  ][1] + E[i-2, j  , k  ][1]))
      + cz.M[0][4] * ((E[i  , j+2, k+1][1] + E[i  , j-2, k+1][1])
                    - (E[i  , j+2, k  ][1] + E[i  , j-2, k  ][1]))
      // n=1: axial k+2, k-1
      + cz.M[1][0] *  (E[i  , j  , k+2][1] - E[i  , j  , k-1][1])
      + cz.M[1][1] * ((E[i+1, j  , k+2][1] + E[i-1, j  , k+2][1])
                    - (E[i+1, j  , k-1][1] + E[i-1, j  , k-1][1]))
      + cz.M[1][2] * ((E[i  , j+1, k+2][1] + E[i  , j-1, k+2][1])
                    - (E[i  , j+1, k-1][1] + E[i  , j-1, k-1][1]))
      + cz.M[1][3] * ((E[i+2, j  , k+2][1] + E[i-2, j  , k+2][1])
                    - (E[i+2, j  , k-1][1] + E[i-2, j  , k-1][1]))
      + cz.M[1][4] * ((E[i  , j+2, k+2][1] + E[i  , j-2, k+2][1])
                    - (E[i  , j+2, k-1][1] + E[i  , j-2, k-1][1]))
      // n=2: axial k+3, k-2
      + cz.M[2][0] *  (E[i  , j  , k+3][1] - E[i  , j  , k-2][1])
      + cz.M[2][1] * ((E[i+1, j  , k+3][1] + E[i-1, j  , k+3][1])
                    - (E[i+1, j  , k-2][1] + E[i-1, j  , k-2][1]))
      + cz.M[2][2] * ((E[i  , j+1, k+3][1] + E[i  , j-1, k+3][1])
                    - (E[i  , j+1, k-2][1] + E[i  , j-1, k-2][1]))
      + cz.M[2][3] * ((E[i+2, j  , k+3][1] + E[i-2, j  , k+3][1])
                    - (E[i+2, j  , k-2][1] + E[i-2, j  , k-2][1]))
      + cz.M[2][4] * ((E[i  , j+2, k+3][1] + E[i  , j-2, k+3][1])
                    - (E[i  , j+2, k-2][1] + E[i  , j-2, k-2][1]));

    // --- D*_y Ez (component 2) using cy coefficients ---
    // y axial, perp1=z (cols 1,3), perp2=x (cols 2,4)
    const auto DyEz =
      // n=0: axial j+1, j
        cy.M[0][0] *  (E[i  , j+1, k  ][2] - E[i  , j  , k  ][2])
      + cy.M[0][1] * ((E[i  , j+1, k+1][2] + E[i  , j+1, k-1][2])
                    - (E[i  , j  , k+1][2] + E[i  , j  , k-1][2]))
      + cy.M[0][2] * ((E[i+1, j+1, k  ][2] + E[i-1, j+1, k  ][2])
                    - (E[i+1, j  , k  ][2] + E[i-1, j  , k  ][2]))
      + cy.M[0][3] * ((E[i  , j+1, k+2][2] + E[i  , j+1, k-2][2])
                    - (E[i  , j  , k+2][2] + E[i  , j  , k-2][2]))
      + cy.M[0][4] * ((E[i+2, j+1, k  ][2] + E[i-2, j+1, k  ][2])
                    - (E[i+2, j  , k  ][2] + E[i-2, j  , k  ][2]))
      // n=1: axial j+2, j-1
      + cy.M[1][0] *  (E[i  , j+2, k  ][2] - E[i  , j-1, k  ][2])
      + cy.M[1][1] * ((E[i  , j+2, k+1][2] + E[i  , j+2, k-1][2])
                    - (E[i  , j-1, k+1][2] + E[i  , j-1, k-1][2]))
      + cy.M[1][2] * ((E[i+1, j+2, k  ][2] + E[i-1, j+2, k  ][2])
                    - (E[i+1, j-1, k  ][2] + E[i-1, j-1, k  ][2]))
      + cy.M[1][3] * ((E[i  , j+2, k+2][2] + E[i  , j+2, k-2][2])
                    - (E[i  , j-1, k+2][2] + E[i  , j-1, k-2][2]))
      + cy.M[1][4] * ((E[i+2, j+2, k  ][2] + E[i-2, j+2, k  ][2])
                    - (E[i+2, j-1, k  ][2] + E[i-2, j-1, k  ][2]))
      // n=2: axial j+3, j-2
      + cy.M[2][0] *  (E[i  , j+3, k  ][2] - E[i  , j-2, k  ][2])
      + cy.M[2][1] * ((E[i  , j+3, k+1][2] + E[i  , j+3, k-1][2])
                    - (E[i  , j-2, k+1][2] + E[i  , j-2, k-1][2]))
      + cy.M[2][2] * ((E[i+1, j+3, k  ][2] + E[i-1, j+3, k  ][2])
                    - (E[i+1, j-2, k  ][2] + E[i-1, j-2, k  ][2]))
      + cy.M[2][3] * ((E[i  , j+3, k+2][2] + E[i  , j+3, k-2][2])
                    - (E[i  , j-2, k+2][2] + E[i  , j-2, k-2][2]))
      + cy.M[2][4] * ((E[i+2, j+3, k  ][2] + E[i-2, j+3, k  ][2])
                    - (E[i+2, j-2, k  ][2] + E[i-2, j-2, k  ][2]));

    Bmds[idx][0] = Bmds[idx][0] + dt * (DzEy - DyEz);
  });


  // === By: By += dt * (D*_x Ez - D*_z Ex) ===
  //
  // D*_x Ez: x is axial, perp1=y (cols 1,3), perp2=z (cols 2,4)
  //   axial offsets in i: n=0 -> (i+1, i), n=1 -> (i+2, i-1), n=2 -> (i+3, i-2)
  //
  // D*_z Ex: z is axial, perp1=x (cols 1,3), perp2=y (cols 2,4)
  //   axial offsets in k: n=0 -> (k+1, k), n=1 -> (k+2, k-1), n=2 -> (k+3, k-2)

  w.for_each_index(Bmds, [=](const auto idx) {
    const auto i = idx[0] + h;
    const auto j = idx[1] + h;
    const auto k = idx[2] + h;
    const auto E = Emds_halo;

    // --- D*_x Ez (component 2) using cx coefficients ---
    // x axial, perp1=y (cols 1,3), perp2=z (cols 2,4)
    const auto DxEz =
      // n=0: axial i+1, i
        cx.M[0][0] *  (E[i+1, j  , k  ][2] - E[i  , j  , k  ][2])
      + cx.M[0][1] * ((E[i+1, j+1, k  ][2] + E[i+1, j-1, k  ][2])
                    - (E[i  , j+1, k  ][2] + E[i  , j-1, k  ][2]))
      + cx.M[0][2] * ((E[i+1, j  , k+1][2] + E[i+1, j  , k-1][2])
                    - (E[i  , j  , k+1][2] + E[i  , j  , k-1][2]))
      + cx.M[0][3] * ((E[i+1, j+2, k  ][2] + E[i+1, j-2, k  ][2])
                    - (E[i  , j+2, k  ][2] + E[i  , j-2, k  ][2]))
      + cx.M[0][4] * ((E[i+1, j  , k+2][2] + E[i+1, j  , k-2][2])
                    - (E[i  , j  , k+2][2] + E[i  , j  , k-2][2]))
      // n=1: axial i+2, i-1
      + cx.M[1][0] *  (E[i+2, j  , k  ][2] - E[i-1, j  , k  ][2])
      + cx.M[1][1] * ((E[i+2, j+1, k  ][2] + E[i+2, j-1, k  ][2])
                    - (E[i-1, j+1, k  ][2] + E[i-1, j-1, k  ][2]))
      + cx.M[1][2] * ((E[i+2, j  , k+1][2] + E[i+2, j  , k-1][2])
                    - (E[i-1, j  , k+1][2] + E[i-1, j  , k-1][2]))
      + cx.M[1][3] * ((E[i+2, j+2, k  ][2] + E[i+2, j-2, k  ][2])
                    - (E[i-1, j+2, k  ][2] + E[i-1, j-2, k  ][2]))
      + cx.M[1][4] * ((E[i+2, j  , k+2][2] + E[i+2, j  , k-2][2])
                    - (E[i-1, j  , k+2][2] + E[i-1, j  , k-2][2]))
      // n=2: axial i+3, i-2
      + cx.M[2][0] *  (E[i+3, j  , k  ][2] - E[i-2, j  , k  ][2])
      + cx.M[2][1] * ((E[i+3, j+1, k  ][2] + E[i+3, j-1, k  ][2])
                    - (E[i-2, j+1, k  ][2] + E[i-2, j-1, k  ][2]))
      + cx.M[2][2] * ((E[i+3, j  , k+1][2] + E[i+3, j  , k-1][2])
                    - (E[i-2, j  , k+1][2] + E[i-2, j  , k-1][2]))
      + cx.M[2][3] * ((E[i+3, j+2, k  ][2] + E[i+3, j-2, k  ][2])
                    - (E[i-2, j+2, k  ][2] + E[i-2, j-2, k  ][2]))
      + cx.M[2][4] * ((E[i+3, j  , k+2][2] + E[i+3, j  , k-2][2])
                    - (E[i-2, j  , k+2][2] + E[i-2, j  , k-2][2]));

    // --- D*_z Ex (component 0) using cz coefficients ---
    // z axial, perp1=x (cols 1,3), perp2=y (cols 2,4)
    const auto DzEx =
      // n=0: axial k+1, k
        cz.M[0][0] *  (E[i  , j  , k+1][0] - E[i  , j  , k  ][0])
      + cz.M[0][1] * ((E[i+1, j  , k+1][0] + E[i-1, j  , k+1][0])
                    - (E[i+1, j  , k  ][0] + E[i-1, j  , k  ][0]))
      + cz.M[0][2] * ((E[i  , j+1, k+1][0] + E[i  , j-1, k+1][0])
                    - (E[i  , j+1, k  ][0] + E[i  , j-1, k  ][0]))
      + cz.M[0][3] * ((E[i+2, j  , k+1][0] + E[i-2, j  , k+1][0])
                    - (E[i+2, j  , k  ][0] + E[i-2, j  , k  ][0]))
      + cz.M[0][4] * ((E[i  , j+2, k+1][0] + E[i  , j-2, k+1][0])
                    - (E[i  , j+2, k  ][0] + E[i  , j-2, k  ][0]))
      // n=1: axial k+2, k-1
      + cz.M[1][0] *  (E[i  , j  , k+2][0] - E[i  , j  , k-1][0])
      + cz.M[1][1] * ((E[i+1, j  , k+2][0] + E[i-1, j  , k+2][0])
                    - (E[i+1, j  , k-1][0] + E[i-1, j  , k-1][0]))
      + cz.M[1][2] * ((E[i  , j+1, k+2][0] + E[i  , j-1, k+2][0])
                    - (E[i  , j+1, k-1][0] + E[i  , j-1, k-1][0]))
      + cz.M[1][3] * ((E[i+2, j  , k+2][0] + E[i-2, j  , k+2][0])
                    - (E[i+2, j  , k-1][0] + E[i-2, j  , k-1][0]))
      + cz.M[1][4] * ((E[i  , j+2, k+2][0] + E[i  , j-2, k+2][0])
                    - (E[i  , j+2, k-1][0] + E[i  , j-2, k-1][0]))
      // n=2: axial k+3, k-2
      + cz.M[2][0] *  (E[i  , j  , k+3][0] - E[i  , j  , k-2][0])
      + cz.M[2][1] * ((E[i+1, j  , k+3][0] + E[i-1, j  , k+3][0])
                    - (E[i+1, j  , k-2][0] + E[i-1, j  , k-2][0]))
      + cz.M[2][2] * ((E[i  , j+1, k+3][0] + E[i  , j-1, k+3][0])
                    - (E[i  , j+1, k-2][0] + E[i  , j-1, k-2][0]))
      + cz.M[2][3] * ((E[i+2, j  , k+3][0] + E[i-2, j  , k+3][0])
                    - (E[i+2, j  , k-2][0] + E[i-2, j  , k-2][0]))
      + cz.M[2][4] * ((E[i  , j+2, k+3][0] + E[i  , j-2, k+3][0])
                    - (E[i  , j+2, k-2][0] + E[i  , j-2, k-2][0]));

    Bmds[idx][1] = Bmds[idx][1] + dt * (DxEz - DzEx);
  });


  // === Bz: Bz += dt * (D*_y Ex - D*_x Ey) ===
  //
  // D*_y Ex: y is axial, perp1=z (cols 1,3), perp2=x (cols 2,4)
  //   axial offsets in j
  //
  // D*_x Ey: x is axial, perp1=y (cols 1,3), perp2=z (cols 2,4)
  //   axial offsets in i

  w.for_each_index(Bmds, [=](const auto idx) {
    const auto i = idx[0] + h;
    const auto j = idx[1] + h;
    const auto k = idx[2] + h;
    const auto E = Emds_halo;

    // --- D*_y Ex (component 0) using cy coefficients ---
    // y axial, perp1=z (cols 1,3), perp2=x (cols 2,4)
    const auto DyEx =
      // n=0: axial j+1, j
        cy.M[0][0] *  (E[i  , j+1, k  ][0] - E[i  , j  , k  ][0])
      + cy.M[0][1] * ((E[i  , j+1, k+1][0] + E[i  , j+1, k-1][0])
                    - (E[i  , j  , k+1][0] + E[i  , j  , k-1][0]))
      + cy.M[0][2] * ((E[i+1, j+1, k  ][0] + E[i-1, j+1, k  ][0])
                    - (E[i+1, j  , k  ][0] + E[i-1, j  , k  ][0]))
      + cy.M[0][3] * ((E[i  , j+1, k+2][0] + E[i  , j+1, k-2][0])
                    - (E[i  , j  , k+2][0] + E[i  , j  , k-2][0]))
      + cy.M[0][4] * ((E[i+2, j+1, k  ][0] + E[i-2, j+1, k  ][0])
                    - (E[i+2, j  , k  ][0] + E[i-2, j  , k  ][0]))
      // n=1: axial j+2, j-1
      + cy.M[1][0] *  (E[i  , j+2, k  ][0] - E[i  , j-1, k  ][0])
      + cy.M[1][1] * ((E[i  , j+2, k+1][0] + E[i  , j+2, k-1][0])
                    - (E[i  , j-1, k+1][0] + E[i  , j-1, k-1][0]))
      + cy.M[1][2] * ((E[i+1, j+2, k  ][0] + E[i-1, j+2, k  ][0])
                    - (E[i+1, j-1, k  ][0] + E[i-1, j-1, k  ][0]))
      + cy.M[1][3] * ((E[i  , j+2, k+2][0] + E[i  , j+2, k-2][0])
                    - (E[i  , j-1, k+2][0] + E[i  , j-1, k-2][0]))
      + cy.M[1][4] * ((E[i+2, j+2, k  ][0] + E[i-2, j+2, k  ][0])
                    - (E[i+2, j-1, k  ][0] + E[i-2, j-1, k  ][0]))
      // n=2: axial j+3, j-2
      + cy.M[2][0] *  (E[i  , j+3, k  ][0] - E[i  , j-2, k  ][0])
      + cy.M[2][1] * ((E[i  , j+3, k+1][0] + E[i  , j+3, k-1][0])
                    - (E[i  , j-2, k+1][0] + E[i  , j-2, k-1][0]))
      + cy.M[2][2] * ((E[i+1, j+3, k  ][0] + E[i-1, j+3, k  ][0])
                    - (E[i+1, j-2, k  ][0] + E[i-1, j-2, k  ][0]))
      + cy.M[2][3] * ((E[i  , j+3, k+2][0] + E[i  , j+3, k-2][0])
                    - (E[i  , j-2, k+2][0] + E[i  , j-2, k-2][0]))
      + cy.M[2][4] * ((E[i+2, j+3, k  ][0] + E[i-2, j+3, k  ][0])
                    - (E[i+2, j-2, k  ][0] + E[i-2, j-2, k  ][0]));

    // --- D*_x Ey (component 1) using cx coefficients ---
    // x axial, perp1=y (cols 1,3), perp2=z (cols 2,4)
    const auto DxEy =
      // n=0: axial i+1, i
        cx.M[0][0] *  (E[i+1, j  , k  ][1] - E[i  , j  , k  ][1])
      + cx.M[0][1] * ((E[i+1, j+1, k  ][1] + E[i+1, j-1, k  ][1])
                    - (E[i  , j+1, k  ][1] + E[i  , j-1, k  ][1]))
      + cx.M[0][2] * ((E[i+1, j  , k+1][1] + E[i+1, j  , k-1][1])
                    - (E[i  , j  , k+1][1] + E[i  , j  , k-1][1]))
      + cx.M[0][3] * ((E[i+1, j+2, k  ][1] + E[i+1, j-2, k  ][1])
                    - (E[i  , j+2, k  ][1] + E[i  , j-2, k  ][1]))
      + cx.M[0][4] * ((E[i+1, j  , k+2][1] + E[i+1, j  , k-2][1])
                    - (E[i  , j  , k+2][1] + E[i  , j  , k-2][1]))
      // n=1: axial i+2, i-1
      + cx.M[1][0] *  (E[i+2, j  , k  ][1] - E[i-1, j  , k  ][1])
      + cx.M[1][1] * ((E[i+2, j+1, k  ][1] + E[i+2, j-1, k  ][1])
                    - (E[i-1, j+1, k  ][1] + E[i-1, j-1, k  ][1]))
      + cx.M[1][2] * ((E[i+2, j  , k+1][1] + E[i+2, j  , k-1][1])
                    - (E[i-1, j  , k+1][1] + E[i-1, j  , k-1][1]))
      + cx.M[1][3] * ((E[i+2, j+2, k  ][1] + E[i+2, j-2, k  ][1])
                    - (E[i-1, j+2, k  ][1] + E[i-1, j-2, k  ][1]))
      + cx.M[1][4] * ((E[i+2, j  , k+2][1] + E[i+2, j  , k-2][1])
                    - (E[i-1, j  , k+2][1] + E[i-1, j  , k-2][1]))
      // n=2: axial i+3, i-2
      + cx.M[2][0] *  (E[i+3, j  , k  ][1] - E[i-2, j  , k  ][1])
      + cx.M[2][1] * ((E[i+3, j+1, k  ][1] + E[i+3, j-1, k  ][1])
                    - (E[i-2, j+1, k  ][1] + E[i-2, j-1, k  ][1]))
      + cx.M[2][2] * ((E[i+3, j  , k+1][1] + E[i+3, j  , k-1][1])
                    - (E[i-2, j  , k+1][1] + E[i-2, j  , k-1][1]))
      + cx.M[2][3] * ((E[i+3, j+2, k  ][1] + E[i+3, j-2, k  ][1])
                    - (E[i-2, j+2, k  ][1] + E[i-2, j-2, k  ][1]))
      + cx.M[2][4] * ((E[i+3, j  , k+2][1] + E[i+3, j  , k-2][1])
                    - (E[i-2, j  , k+2][1] + E[i-2, j  , k-2][1]));

    Bmds[idx][2] = Bmds[idx][2] + dt * (DyEx - DxEy);
  });
}

}  // namespace emf
