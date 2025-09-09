#include "core/pic2/particle.h"

#include "core/mdgrid_common.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/sort.h"
#include "thrust/unique.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <algorithm>
#include <cmath>
#include <iterator>
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

void
  ParticleContainer::wrap_positions(
    const std::array<value_type, 3> mins,
    const std::array<value_type, 3> maxs)
{
  if(mins[0] >= maxs[0] or mins[1] >= maxs[1] or mins[0] >= maxs[2]) {
    throw std::runtime_error {
      "ParticleContainer::wrap_coordinates: given limits do not make sense."
    };
  }
  const auto L = std::array { maxs[0] - mins[0], maxs[1] - mins[1], maxs[2] - mins[2] };

  const auto pos_mds = this->pos_.mds();
  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx, const auto tidx) {
        const auto a = mins[tidx[0]];
        const auto b = maxs[tidx[0]];
        const auto x = pos_mds[idx][tidx];

        const auto wrapped_pos = (x < 0 ? b : a) + std::fmod(x, L[tidx[0]]);
        pos_mds[idx][tidx]     = wrapped_pos;
      })
    .wait();
}

std::pair<std::map<std::array<int, 3>, ParticleContainer::span>, ParticleContainer>
  ParticleContainer::divide_to_subregions(
    const std::array<value_type, 2> x_dividers,
    const std::array<value_type, 2> y_dividers,
    const std::array<value_type, 2> z_dividers)
{
  if(
    x_dividers[0] >= x_dividers[1] or y_dividers[0] >= y_dividers[1] or
    z_dividers[0] >= z_dividers[1]) {
    throw std::runtime_error {
      "ParticleContainer::split_to_subregions: given dividers do not make sense."
    };
  }

  // Handle empty container as special case, so that we can later assume
  // that there is at least one particle.
  if(this->size() == 0) { return { {}, *this }; }

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

  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();

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

  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<std::size_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->size());

  const auto particle_subregion_indices_begin = thrust::make_transform_iterator(
    particle_ordinals_begin,
    particle_ordinal_to_subregion_index);
  const auto particle_subregion_indices_end =
    rn::next(particle_subregion_indices_begin, this->size());

  auto particle_subregion_indices = thrust::device_vector<std::size_t>(
    particle_subregion_indices_begin,
    particle_subregion_indices_end);

  auto particle_trackers =
    thrust::device_vector<std::size_t>(particle_ordinals_begin, particle_ordinals_end);

  thrust::sort_by_key(
    particle_subregion_indices.begin(),
    particle_subregion_indices.end(),
    particle_trackers.begin());

  auto permuted_pcontainer    = *this;
  const auto permuted_pos_mds = permuted_pcontainer.pos_.mds();
  const auto permuted_vel_mds = permuted_pcontainer.vel_.mds();

  tyvi::mdgrid_work {}
    .for_each_index(
      permuted_pos_mds,
      [=, p = particle_trackers.begin()](const auto idx, const auto tidx) {
        const std::size_t i         = p[idx[0]];
        permuted_pos_mds[idx][tidx] = pos_mds[i][tidx];
        permuted_vel_mds[idx][tidx] = vel_mds[i][tidx];
      })
    .wait();

  // Reset these back to 0, 1, ... inorder to use them in unique
  // to measure the begins of the subregions.
  particle_trackers.assign(particle_ordinals_begin, particle_ordinals_end);

  const auto new_end = thrust::unique_by_key(
    thrust::device,
    particle_subregion_indices.begin(),
    particle_subregion_indices.end(),
    particle_trackers.begin());

  // Note that this ParticleContainer might not contain
  // particles in all of the 27 subregions.

  const auto present_subregion_indices =
    thrust::host_vector<std::size_t>(particle_subregion_indices.begin(), new_end.first);
  const auto present_subregion_begins =
    thrust::host_vector<std::size_t>(particle_trackers.begin(), new_end.second);

  const auto present_subregion_ends = [&] {
    // We assumed at beginning of the function that there is at least one particle.
    auto temp =
      std::vector(++present_subregion_begins.begin(), present_subregion_begins.end());
    temp.push_back(this->size());
    return temp;
  }();

  std::map<std::array<int, 3>, ParticleContainer::span> m {};

  for(auto i = 0uz; i < present_subregion_indices.size(); ++i) {
    m.insert_or_assign(
      subregion_index_to_dir(present_subregion_indices[i]),
      ParticleContainer::span { present_subregion_begins[i],
                                present_subregion_ends[i] });
  }
  return { std::move(m), std::move(permuted_pcontainer) };
}

}  // namespace pic2
