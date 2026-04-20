#include "core/pic/particle.h"

#include "core/mdgrid_common.h"
#include "thrust/count.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/reduce.h"
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

namespace pic {

ParticleContainer::ParticleContainer(const ParticleContainerArgs args) :
  pos_(args.N),
  vel_(args.N),
  ids_(args.N),
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
  auto w = tyvi::mdgrid_work {};
  w.sync_to_staging(pos_);
  w.sync_to_staging(ids_);

  auto x = std::vector<value_type>(this->size());
  auto y = std::vector<value_type>(this->size());
  auto z = std::vector<value_type>(this->size());

  w.wait();

  const auto mds   = pos_.staging_mds();
  const auto idmds = ids_.staging_mds();

  auto n = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    x[n] = mds[i][0];
    y[n] = mds[i][1];
    z[n] = mds[i][2];
    ++n;
  }

  x.resize(n);
  y.resize(n);
  z.resize(n);

  return { std::move(x), std::move(y), std::move(z) };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_velocities()
{

  auto w = tyvi::mdgrid_work {};
  w.sync_to_staging(vel_);
  w.sync_to_staging(ids_);

  auto vx = std::vector<value_type>(this->size());
  auto vy = std::vector<value_type>(this->size());
  auto vz = std::vector<value_type>(this->size());

  w.wait();

  const auto mds   = vel_.staging_mds();
  const auto idmds = ids_.staging_mds();
  auto n           = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    vx[n] = mds[i][0];
    vy[n] = mds[i][1];
    vz[n] = mds[i][2];
    ++n;
  }

  vx.resize(n);
  vy.resize(n);
  vz.resize(n);

  return { std::move(vx), std::move(vy), std::move(vz) };
}

std::vector<runko::prtc_id_type>
  ParticleContainer::get_ids()
{

  auto w = tyvi::mdgrid_work {};
  w.sync_to_staging(ids_);

  auto tmp = std::vector<runko::prtc_id_type>(this->size());

  w.wait();

  const auto idmds = ids_.staging_mds();
  auto n           = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(idmds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    tmp[n] = idmds[i][];
    ++n;
  }

  tmp.resize(n);
  return tmp;
}


void
  ParticleContainer::wrap_positions(
    const std::array<value_type, 3> mins,
    const std::array<value_type, 3> maxs)
{
  if(mins[0] >= maxs[0] or mins[1] >= maxs[1] or mins[2] >= maxs[2]) {
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

std::map<runko::grid_neighbor<3>, std::pair<std::size_t, std::size_t>>
  ParticleContainer::divide_to_subregions(
    thrust::device_vector<runko::ParticleState<ParticleContainer::value_type>>& buffer,
    std::array<value_type, 2> x_dividers,
    std::array<value_type, 2> y_dividers,
    std::array<value_type, 2> z_dividers)
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
    [](const int i, const int j, const int k) -> runko::index_t {
    return runko::grid_neighbor<3>(i, j, k).neighbor_index();
  };

  static constexpr auto subregion_index_to_dir = [](const runko::index_t i) static {
    return runko::grid_neighbor<3>::from_index(i);
  };

  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();
  const auto ids_mds = this->ids_.mds();

  const auto coords_to_subregion_index = [=](const auto x, const auto y, const auto z) {
    const auto i =
      static_cast<int>(x >= x_dividers[0]) - static_cast<int>(x < x_dividers[1]);
    const auto j =
      static_cast<int>(y >= y_dividers[0]) - static_cast<int>(y < y_dividers[1]);
    const auto k =
      static_cast<int>(z >= z_dividers[0]) - static_cast<int>(z < z_dividers[1]);

    return subregion_index(i, j, k);
  };

  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->size());

  const auto is_leaving = [=](const auto index) {
    return index != subregion_index(0, 0, 0);
  };

  tyvi::mdgrid_work w {};
  const auto leaving_particles = thrust::count_if(
    w.on_this(),
    particle_ordinals_begin,
    particle_ordinals_end,
    [=](const auto i) {
      return ids_mds[i][] != runko::dead_prtc_id and
             is_leaving(
               coords_to_subregion_index(pos_mds[i][0], pos_mds[i][1], pos_mds[i][2]));
    });

  const auto prev_buffer_size = buffer.size();
  buffer.resize(prev_buffer_size + leaving_particles);

  const auto particle_states_begin =
    thrust::make_transform_iterator(particle_ordinals_begin, [=](const auto n) {
      return runko::ParticleState<value_type> {
        .pos { pos_mds[n][0], pos_mds[n][1], pos_mds[n][2] },
        .vel { vel_mds[n][0], vel_mds[n][1], vel_mds[n][2] },
        .id { ids_mds[n][] }
      };
    });
  const auto particle_states_end = rn::next(particle_states_begin, this->size());

  const auto buffer_data_begin = rn::next(buffer.begin(), prev_buffer_size);

  const auto copy_if_result = thrust::copy_if(
    w.on_this(),
    particle_states_begin,
    particle_states_end,
    buffer_data_begin,
    [=](const auto& state) {
      return state.id != runko::dead_prtc_id and
             is_leaving(
               coords_to_subregion_index(state.pos[0], state.pos[1], state.pos[2]));
    });

  if(copy_if_result != buffer.end()) {
    throw std::logic_error {
      "These should match, because buffer is resized specificly for correct amount of "
      "particles."
    };
  }

  tyvi::mdgrid_work w2 {};
  tyvi::when_all(w, w2);
  w2.for_each_index(pos_mds, [=](const auto idx) {
    if(is_leaving(coords_to_subregion_index(
         pos_mds[idx][0],
         pos_mds[idx][1],
         pos_mds[idx][2]))) {
      ids_mds[idx][] = runko::dead_prtc_id;
    }
  });

  thrust::sort(
    w.on_this(),
    buffer_data_begin,
    buffer.end(),
    [=](const auto& lhs, const auto& rhs) {
      const auto sublhs = coords_to_subregion_index(lhs.pos[0], lhs.pos[1], lhs.pos[2]);
      const auto subrhs = coords_to_subregion_index(rhs.pos[0], rhs.pos[1], rhs.pos[2]);
      return sublhs < subrhs;
    });

  using subregion_counter     = toolbox::VecD<runko::index_t, 27>;
  const auto state_to_counter = [=](const runko::ParticleState<value_type>& x) {
    const auto i = coords_to_subregion_index(x.pos[0], x.pos[1], x.pos[2]);
    auto v       = subregion_counter {};
    v[i]         = 1;
    return v;
  };

  const auto counters_begin =
    thrust::transform_iterator(buffer_data_begin, state_to_counter);
  const auto counters_end = rn::next(counters_begin, leaving_particles);

  const auto subregion_counts =
    thrust::reduce(w.on_this(), counters_begin, counters_end, subregion_counter {});

  std::map<runko::grid_neighbor<3>, std::pair<std::size_t, std::size_t>> m {};

  auto next_span_begin = prev_buffer_size;
  for(auto i = runko::index_t { 0 }; i < runko::index_t { 27 }; ++i) {
    const auto dir = subregion_index_to_dir(i);
    if(i == subregion_index(0, 0, 0)) {
      m.insert_or_assign(dir, std::pair { next_span_begin, next_span_begin });
      continue;
    }

    const auto count = subregion_counts[i];
    const auto from  = next_span_begin;
    const auto until = from + count;
    next_span_begin  = until;

    m.insert_or_assign(dir, std::pair { from, until });
  }

  tyvi::when_all(w, w2);
  w.wait();
  return m;
}


double
  ParticleContainer::total_kinetic_energy() const
{
  auto w = tyvi::mdgrid_work {};

  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);

  // Note that these are in natural units.
  const auto vel_mds = this->vel_.mds();
  const auto ids_mds = this->ids_.mds();

  auto particle_ordinal_to_kinetic_energy = [=](const runko::index_t n) -> double {
    const auto is_alive = ids_mds[n][] != runko::dead_prtc_id;
    using Vec3          = toolbox::Vec3<value_type>;
    const auto v        = Vec3(vel_mds[n]);
    const auto e = std::sqrt(value_type { 1 } + toolbox::dot(v, v)) - value_type { 1 };
    return is_alive * static_cast<double>(e);
  };

  const auto kinetic_energies_begin = thrust::make_transform_iterator(
    particle_ordinals_begin,
    particle_ordinal_to_kinetic_energy);
  const auto kinetic_energies_end = rn::next(kinetic_energies_begin, this->size());

  return thrust::reduce(w.on_this(), kinetic_energies_begin, kinetic_energies_end);
}

}  // namespace pic
