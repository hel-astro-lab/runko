#include "core/emf/yee_lattice.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <tuple>
#include <utility>

void
  emf::YeeLattice::filter_current_binomial2()
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
    std::tuple { std::integral_constant<std::size_t, 1uz> {}, Nx - 1uz },
    std::tuple { std::integral_constant<std::size_t, 1uz> {}, Ny - 1uz },
    std::tuple { std::integral_constant<std::size_t, 1uz> {}, Nz - 1uz });


  const auto Jmds = this->J_.mds();

  tyvi::mdgrid_work {}
    .for_each_index(
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
