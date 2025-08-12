#include "core/pic2/particle.h"

#include "core/mdgrid_common.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/transform.h"
#include "thrust/transform_scan.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <cmath>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace pic2 {

ParticleContainer::ParticleContainer(const ParticleContainerArgs args) :
  pos_(args.N),
  vel_(args.N),
  charge_ { args.charge },
  mass_ { args.mass }
{
}

std::size_t
  ParticleContainer::size() const
{
  return pos_.extents().extent(0);
}

ParticleContainerArgs
  ParticleContainer::args() const
{
  return { .N = size(), .charge = charge_, .mass = mass_ };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_positions()
{
  auto w = tyvi::mdgrid_work {}.sync_to_staging(pos_);

  auto x = std::vector<value_type>(this->size());
  auto y = std::vector<value_type>(this->size());
  auto z = std::vector<value_type>(this->size());

  w.wait();

  const auto mds = pos_.staging_mds();

  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    x[i] = mds[i][0];
    y[i] = mds[i][1];
    z[i] = mds[i][2];
  }

  return { std::move(x), std::move(y), std::move(z) };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_velocities()
{

  auto w = tyvi::mdgrid_work {}.sync_to_staging(vel_);

  auto vx = std::vector<value_type>(this->size());
  auto vy = std::vector<value_type>(this->size());
  auto vz = std::vector<value_type>(this->size());

  w.wait();

  const auto mds = vel_.staging_mds();

  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    vx[i] = mds[i][0];
    vy[i] = mds[i][1];
    vz[i] = mds[i][2];
  }

  return { std::move(vx), std::move(vy), std::move(vz) };
}

std::vector<std::pair<std::array<int, 3>, ParticleContainer>>
  ParticleContainer::split_to_subregions(
    const std::array<value_type, 2> x_dividers,
    const std::array<value_type, 2> y_dividers,
    const std::array<value_type, 2> z_dividers,
    const std::array<value_type, 3> global_mins,
    const std::array<value_type, 3> global_maxs)
{
  if(
    x_dividers[0] >= x_dividers[1] or y_dividers[0] >= y_dividers[1] or
    z_dividers[0] >= z_dividers[1]) {
    throw std::runtime_error {
      "ParticleContainer::split_to_subregions: given dividers do not make sense."
    };
  }

  // subregion_index: {-1, 0, 1} x {-1, 0, 1} x {-1, 0, 1} -> {0, 1, ..., 26}
  constexpr auto subregion_index =
    [](const int i, const int j, const int k) -> std::size_t {
    const auto ii = static_cast<std::size_t>(i + 1);
    const auto jj = static_cast<std::size_t>(j + 1);
    const auto kk = static_cast<std::size_t>(k + 1);
    return ii * 9 + jj * 3 + kk;
  };

  constexpr auto subregion_index_to_dir = [](const std::size_t nuz) {
    const auto n = static_cast<int>(nuz);
    return std::array<int, 3> { (n / 9) - 1, ((n / 3) % 3) - 1, (n % 3) - 1 };
  };

  // Special case of zero particles.
  if(this->size() == 0) {
    auto empty_subregion_pcontainers =
      std::vector<std::pair<std::array<int, 3>, ParticleContainer>> {};
    empty_subregion_pcontainers.reserve(27uz);

    const auto args             = this->args();
    const auto empty_pcontainer = ParticleContainer(args);


    for(const auto i: std::views::iota(0uz, 27uz)) {
      const auto dir = subregion_index_to_dir(i);
      empty_subregion_pcontainers.push_back({ dir, empty_pcontainer });
    }

    static constexpr auto my_idx = subregion_index(0, 0, 0);
    empty_subregion_pcontainers.erase(empty_subregion_pcontainers.begin() + my_idx);

    return empty_subregion_pcontainers;
  }

  const auto pos_mds = this->pos_.mds();
  const auto particle_ordinal_to_subregion_index =
    [=](const std::size_t n) -> std::size_t {
    const auto x = pos_mds[n][0];
    const auto y = pos_mds[n][1];
    const auto z = pos_mds[n][2];

    const auto i =
      static_cast<int>(x >= x_dividers[0]) - static_cast<int>(x < x_dividers[1]);
    const auto j =
      static_cast<int>(y >= y_dividers[0]) - static_cast<int>(y < y_dividers[1]);
    const auto k =
      static_cast<int>(z >= z_dividers[0]) - static_cast<int>(z < z_dividers[1]);

    return subregion_index(i, j, k);
  };

  const auto particle_ordinals_begin = thrust::counting_iterator<std::size_t>(0uz);
  const auto particle_ordinals_end   = particle_ordinals_begin + this->size();
  auto subregion_indices             = thrust::device_vector<std::size_t>(this->size());

  thrust::transform(
    thrust::device,
    particle_ordinals_begin,
    particle_ordinals_end,
    subregion_indices.begin(),
    particle_ordinal_to_subregion_index);

  using SubRegionSizes = toolbox::VecD<std::size_t, 27>;

  const auto subregion_index_to_size_contribution =
    [=](const std::size_t i) -> SubRegionSizes {
    auto c = SubRegionSizes {};
    c(i)   = 1uz;
    return c;
  };

  auto subregion_sizes = thrust::device_vector<SubRegionSizes>(this->size());

  thrust::transform_inclusive_scan(
    thrust::device,
    subregion_indices.begin(),
    subregion_indices.end(),
    subregion_sizes.begin(),
    subregion_index_to_size_contribution,
    thrust::plus<SubRegionSizes> {});

  // Special case of zero particles is handled in the beginning.
  const SubRegionSizes final_sizes = subregion_sizes.back();

  auto subregion_pcontainers =
    std::vector<std::pair<std::array<int, 3>, ParticleContainer>> {};
  subregion_pcontainers.reserve(27uz);

  for(const auto i: std::views::iota(0uz, 27uz)) {
    auto a = this->args();
    a.N    = final_sizes(i);

    const auto dir = subregion_index_to_dir(i);
    subregion_pcontainers.push_back({ dir, ParticleContainer(a) });
  }

  using VecMDS =
    std::remove_cvref_t<decltype(std::declval<decltype(this->pos_)&>().mds())>;

  auto subregion_pos_mds = std::array<VecMDS, 27> {};
  auto subregion_vel_mds = std::array<VecMDS, 27> {};

  for(const auto i: std::views::iota(0uz, 27uz)) {
    subregion_pos_mds[i] = subregion_pcontainers[i].second.pos_.mds();
    subregion_vel_mds[i] = subregion_pcontainers[i].second.vel_.mds();
  }

  const auto vel_mds               = this->vel_.mds();
  const auto subregion_indices_ptr = subregion_indices.begin();
  const auto subregion_sizes_ptr   = subregion_sizes.begin();

  const auto global_L = std::array { global_maxs[0] - global_mins[0],
                                     global_maxs[1] - global_mins[1],
                                     global_maxs[2] - global_mins[2] };

  tyvi::mdgrid_work {}
    .for_each_index(
      vel_mds,
      [=](const auto idx, const auto tidx) {
        const auto subregion_index            = subregion_indices_ptr[idx[0]];
        const SubRegionSizes& subregion_sizes = subregion_sizes_ptr[idx[0]];
        const auto offset_in_subregion        = subregion_sizes(subregion_index) - 1uz;

        const auto a = global_mins[tidx[0]];
        const auto b = global_maxs[tidx[0]];
        const auto x = pos_mds[idx][tidx];

        const auto wrapped_pos = (x < 0 ? b : a) + std::fmod(x, global_L[tidx[0]]);

        subregion_pos_mds[subregion_index][offset_in_subregion][tidx] = wrapped_pos;
        subregion_vel_mds[subregion_index][offset_in_subregion][tidx] =
          vel_mds[idx][tidx];
      })
    .wait();

  static constexpr auto my_idx = subregion_index(0, 0, 0);
  *this                        = std::move(subregion_pcontainers[my_idx].second);

  subregion_pcontainers.erase(subregion_pcontainers.begin() + my_idx);

  return subregion_pcontainers;
}

}  // namespace pic2
