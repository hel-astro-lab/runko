#include "core/emf/yee_lattice.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <tuple>
#include <utility>

void
  emf::YeeLattice::filter_current_binomial2()
{
  const tyvi::mdgrid_work w {};
  filter_current_binomial2(w);
}

// 3D single-pass binomial filter using a 3x3x3 coefficient array.
// NOTE: The C3 array is found to block CPU auto-SIMD vectorization;
//       prefer the separable binomial2 variant for CPU backends.

void
  emf::YeeLattice::filter_current_binomial2_3d()
{
  const tyvi::mdgrid_work w {};
  filter_current_binomial2_3d(w);
}

void
  emf::YeeLattice::filter_current_binomial2_3d(const tyvi::mdgrid_work& w)
{
  // 3D 3-point binomial coefficients
  static constexpr value_type C3[3][3][3] = { { { 1. / 64., 2. / 64., 1. / 64. },
                                                { 2. / 64., 4. / 64., 2. / 64. },
                                                { 1. / 64., 2. / 64., 1. / 64. } },
                                              { { 2. / 64., 4. / 64., 2. / 64. },
                                                { 4. / 64., 8. / 64., 4. / 64. },
                                                { 2. / 64., 4. / 64., 2. / 64. } },
                                              { { 1. / 64., 2. / 64., 1. / 64. },
                                                { 2. / 64., 4. / 64., 2. / 64. },
                                                { 1. / 64., 2. / 64., 1. / 64. } } };
  constexpr auto C3_index_space           = tyvi::sstd::geometric_index_space<3, 3>();

  auto filteredJ           = VecGrid(this->J_.extents());
  const auto e             = filteredJ.extents();
  const auto Nx            = e.extent(0);
  const auto Ny            = e.extent(1);
  const auto Nz            = e.extent(2);
  const auto filteredJ_mds = std::submdspan(
    filteredJ.mds(),
    std::tuple { std::integral_constant<runko::index_t, 1u> {},
                 static_cast<runko::index_t>(Nx - 1u) },
    std::tuple { std::integral_constant<runko::index_t, 1u> {},
                 static_cast<runko::index_t>(Ny - 1u) },
    std::tuple { std::integral_constant<runko::index_t, 1u> {},
                 static_cast<runko::index_t>(Nz - 1u) });


  const auto Jmds = this->J_.mds();

  w.for_each_index(
     filteredJ_mds,
     [=](const auto idx, const auto tidx) {
       filteredJ_mds[idx][tidx] = 0;

       for(const auto jdx: C3_index_space) {
         const auto i = idx[0] + jdx[0];
         const auto j = idx[1] + jdx[1];
         const auto k = idx[2] + jdx[2];

         // operator += is broken in hip
         filteredJ_mds[idx][tidx] =
           filteredJ_mds[idx][tidx] + C3[jdx[0]][jdx[1]][jdx[2]] * Jmds[i, j, k][tidx];
       }
     })
    .wait();

  this->J_ = std::move(filteredJ);
}


void
  emf::YeeLattice::filter_current_binomial2(const tyvi::mdgrid_work& w)
{
  // 1D binomial coefficients (the 3D stencil is separable: C3[a][b][c] =
  // B1[a]*B1[b]*B1[c])
  static constexpr value_type B1[3] = { 0.25f, 0.5f, 0.25f };

  const auto e  = this->J_.extents();
  const auto Nx = e.extent(0);
  const auto Ny = e.extent(1);
  const auto Nz = e.extent(2);

  const auto Jmds = this->J_.mds();

  // Pass 1: convolve along z (dim 2)
  // temp1 dims: (Nx, Ny, Nz-2) — z-dimension shrinks by 2
  auto temp1    = VecGrid(Nx, Ny, Nz - 2uz);
  const auto t1 = temp1.mds();

  w.for_each_index(
     t1,
     [=](const auto idx, const auto tidx) {
       const auto i = idx[0], j = idx[1], k = idx[2];
       t1[idx][tidx] = B1[0] * Jmds[i, j, k][tidx] + B1[1] * Jmds[i, j, k + 1][tidx] +
                       B1[2] * Jmds[i, j, k + 2][tidx];
     })
    .wait();

  // Pass 2: convolve along y (dim 1)
  // temp2 dims: (Nx, Ny-2, Nz-2) — y-dimension shrinks by 2
  auto temp2    = VecGrid(Nx, Ny - 2uz, Nz - 2uz);
  const auto t2 = temp2.mds();

  w.for_each_index(
     t2,
     [=](const auto idx, const auto tidx) {
       const auto i = idx[0], j = idx[1], k = idx[2];
       t2[idx][tidx] = B1[0] * t1[i, j, k][tidx] + B1[1] * t1[i, j + 1, k][tidx] +
                       B1[2] * t1[i, j + 2, k][tidx];
     })
    .wait();

  // Pass 3: convolve along x (dim 0) — write to interior of full-size output
  auto filteredJ           = VecGrid(this->J_.extents());
  const auto filteredJ_mds = std::submdspan(
    filteredJ.mds(),
    std::tuple { std::integral_constant<runko::index_t, 1uz> {},
                 static_cast<runko::index_t>(Nx - 1) },
    std::tuple { std::integral_constant<runko::index_t, 1uz> {},
                 static_cast<runko::index_t>(Ny - 1) },
    std::tuple { std::integral_constant<runko::index_t, 1uz> {},
                 static_cast<runko::index_t>(Nz - 1) });

  w.for_each_index(
     filteredJ_mds,
     [=](const auto idx, const auto tidx) {
       const auto i = idx[0], j = idx[1], k = idx[2];
       filteredJ_mds[idx][tidx] = B1[0] * t2[i, j, k][tidx] +
                                  B1[1] * t2[i + 1, j, k][tidx] +
                                  B1[2] * t2[i + 2, j, k][tidx];
     })
    .wait();

  this->J_ = std::move(filteredJ);
}
