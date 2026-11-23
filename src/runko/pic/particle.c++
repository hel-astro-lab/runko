// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/pic/particle.h"

#include "runko/mdgrid_common.h"
#include "runko/tools/vector.h"
#include "thrust/copy.h"
#include "thrust/count.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/fill.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/reduce.h"
#include "thrust/sort.h"
#include "thrust/transform.h"
#include "thrust/unique.h"
#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace pic {

ParticleContainer::ParticleContainer(const ParticleContainerArgs args) :
  pos_(args.N),
  vel_(args.N),
  ids_(args.N),
  wgt_(args.mass == 0.0 ? args.N : 0uz),
  payload_(args.n_payload),
  charge_ { args.charge },
  mass_ { args.mass },
  weight_ { args.weight }
{
  for(auto& p: payload_) { p.resize(args.N); }
}

std::size_t
ParticleContainer::size() const
{
  return pos_.size();
}

std::ptrdiff_t
  ParticleContainer::ssize() const
{
  return static_cast<std::ptrdiff_t>(this->size());
}

std::size_t
  ParticleContainer::number_alive() const
{
  const auto ids_mds = this->ids_.mds();
  const auto i0      = thrust::counting_iterator<runko::index_t>(0);
  const auto w       = tyvi::mdgrid_work {};
  return static_cast<std::size_t>(thrust::count_if(
    w.on_this(), i0, i0 + this->ssize(), [=](const runko::index_t i) {
      return ids_mds[i][] != runko::dead_prtc_id;
    }));
}


ParticleContainerArgs
  ParticleContainer::args() const
{
  return { .N                    = size(),
           .charge               = charge_,
           .mass                 = mass_,
           .weight               = weight_,
           .n_payload            = n_payload() };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_positions()
{
  auto host_pos = runko::host_vec_segments<value_type>(this->size());
  auto host_ids = runko::host_scalar_segments<runko::prtc_id_type>(this->size());

  tyvi::copy_d2h(this->pos_, host_pos);
  tyvi::copy_d2h(this->ids_, host_ids);

  auto x = std::vector<value_type>(this->size());
  auto y = std::vector<value_type>(this->size());
  auto z = std::vector<value_type>(this->size());

  const auto mds   = host_pos.mds();
  const auto idmds = host_ids.mds();

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

  auto host_vel = runko::host_vec_segments<value_type>(this->size());
  auto host_ids = runko::host_scalar_segments<runko::prtc_id_type>(this->size());

  tyvi::copy_d2h(this->vel_, host_vel);
  tyvi::copy_d2h(this->ids_, host_ids);

  auto vx = std::vector<value_type>(this->size());
  auto vy = std::vector<value_type>(this->size());
  auto vz = std::vector<value_type>(this->size());

  const auto mds   = host_vel.mds();
  const auto idmds = host_ids.mds();

  auto n = 0uz;
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

  auto host_ids = runko::host_scalar_segments<runko::prtc_id_type>(this->size());
  tyvi::copy_d2h(this->ids_, host_ids);

  auto tmp = std::vector<runko::prtc_id_type>(this->size());

  const auto idmds = host_ids.mds();
  auto n           = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(idmds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    tmp[n] = idmds[i][];
    ++n;
  }

  tmp.resize(n);
  return tmp;
}

std::vector<ParticleContainer::value_type>
  ParticleContainer::get_weights()
{
  auto host_ids = runko::host_scalar_segments<runko::prtc_id_type>(this->size());
  tyvi::copy_d2h(this->ids_, host_ids);
  auto host_wgt = runko::host_scalar_segments<value_type>(this->wgt_.size());
  if(this->has_particle_weights()) { tyvi::copy_d2h(this->wgt_, host_wgt); }

  auto tmp = std::vector<value_type>(this->size());

  const auto idmds = host_ids.mds();
  const auto wmds  = host_wgt.mds();
  auto n           = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(idmds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    tmp[n] = this->has_particle_weights() ? wmds[i][] : static_cast<value_type>(weight_);
    ++n;
  }

  tmp.resize(n);
  return tmp;
}

std::vector<ParticleContainer::value_type>
  ParticleContainer::get_payload(const std::size_t c)
{
  auto host_ids = runko::host_scalar_segments<runko::prtc_id_type>(this->size());
  auto host_p   = runko::host_scalar_segments<value_type>(this->size());
  tyvi::copy_d2h(this->ids_, host_ids);
  tyvi::copy_d2h(this->payload_.at(c), host_p);

  auto tmp = std::vector<value_type>(this->size());

  const auto idmds = host_ids.mds();
  const auto pmds  = host_p.mds();
  auto n           = 0uz;
  for(const auto [i]: tyvi::sstd::index_space(idmds)) {
    if(idmds[i][] == runko::dead_prtc_id) { continue; }

    tmp[n] = pmds[i][];
    ++n;
  }

  tmp.resize(n);
  return tmp;
}

void
  ParticleContainer::reorder(const thrust::device_vector<runko::index_t>& src)
{
  /* The component views are spelled out because of a bug in gcc:
     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=12526 */
  using I                    = runko::index_t;
  static constexpr auto zero = std::array { I { 0 } };
  static constexpr auto one  = std::array { I { 1 } };
  static constexpr auto two  = std::array { I { 2 } };

  const auto n_old      = this->size();
  const auto n_new      = src.size();
  const auto resize_all = [&](const std::size_t n) {
    this->pos_.resize(n);
    this->vel_.resize(n);
    this->ids_.resize(n);
    if(this->has_particle_weights()) { this->wgt_.resize(n); }
    for(auto& p: this->payload_) { p.resize(n); }
  };
  if(n_new > n_old) { resize_all(n_new); }

  // Sized once for the widest component, so no reallocation while work is in flight.
  this->generic_cache_.resize(
    std::max(sizeof(value_type), sizeof(runko::prtc_id_type)) * n_old);

  const auto w = tyvi::mdgrid_work {};
  // Copy the first n_old values of one component aside, then write slot j from old slot
  // src[j]. The components share the cache, so their work is sequenced on w.
  const auto move = [&]<typename T>(const auto comp, const T fill) {
    const auto p =
      reinterpret_cast<T*>(thrust::raw_pointer_cast(this->generic_cache_.data()));
    if(reinterpret_cast<std::uintptr_t>(p) % alignof(T)) {
      throw std::runtime_error { "ParticleContainer::reorder: misaligned cache." };
    }
    const auto first = comp.begin();
    std::ignore      = thrust::copy(
      w.on_this(),
      first,
      std::ranges::next(first, static_cast<std::ptrdiff_t>(n_old)),
      p);
    std::ignore = thrust::transform(
      w.on_this(),
      src.begin(),
      src.end(),
      first,
      [=](const I s) { return s == no_source ? fill : p[s]; });
  };

  // pos/vel/wgt/payload of a dead slot are never read: only the id distinguishes it.
  move(this->pos_.component_view<zero>(), value_type { 0 });
  move(this->pos_.component_view<one>(), value_type { 0 });
  move(this->pos_.component_view<two>(), value_type { 0 });
  move(this->vel_.component_view<zero>(), value_type { 0 });
  move(this->vel_.component_view<one>(), value_type { 0 });
  move(this->vel_.component_view<two>(), value_type { 0 });
  if(this->has_particle_weights()) {
    move(this->wgt_.component_view<>(), value_type { 0 });
  }
  for(auto& p: this->payload_) { move(p.component_view<>(), value_type { 0 }); }
  move(this->ids_.component_view<>(), runko::dead_prtc_id);
  w.wait();

  if(n_new < n_old) { resize_all(n_new); }
  ++this->reorder_count_;
}

void
  ParticleContainer::append_dead(const std::size_t n)
{
  if(n == 0uz) { return; }
  const auto old = this->size();
  this->pos_.resize(old + n);
  this->vel_.resize(old + n);
  this->ids_.resize(old + n);
  if(this->has_particle_weights()) { this->wgt_.resize(old + n); }
  for(auto& p: this->payload_) { p.resize(old + n); }

  // Only the id distinguishes a dead slot. Pushers, depositors, and reflectors
  // early-return on id == dead_prtc_id, so pos_ and vel_ are never read for
  // dead particles and are left uninitialized.
  const auto w        = tyvi::mdgrid_work {};
  const auto ids_span = this->ids_.component_view<>();
  thrust::fill(
    w.on_this(),
    std::ranges::next(ids_span.begin(), static_cast<std::ptrdiff_t>(old)),
    ids_span.end(),
    runko::dead_prtc_id);
  w.wait();
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
    thrust::device_vector<value_type>& payload_buffer,
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

  const auto pos_mds  = this->pos_.mds();
  const auto vel_mds  = this->vel_.mds();
  const auto ids_mds  = this->ids_.mds();
  const auto wgt_mds  = this->wgt_.mds();
  const auto weighted = this->has_particle_weights();
  const auto w0       = static_cast<value_type>(this->weight_);

  // subregion of the particle in slot n
  const auto slot_to_subregion_index = [=](const runko::index_t n) {
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
  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->ssize());

  // The slots of the leaving particles, ordered by subregion and within one by slot.
  tyvi::mdgrid_work w {};
  auto& leaving_slots = this->tracker_cache_;
  leaving_slots.resize(this->size());
  const auto leaving_end = thrust::copy_if(
    w.on_this(),
    particle_ordinals_begin,
    particle_ordinals_end,
    leaving_slots.begin(),
    [=](const auto n) {
      return ids_mds[n][] != runko::dead_prtc_id and
             slot_to_subregion_index(n) != subregion_index(0, 0, 0);
    });
  const auto leaving_particles = rn::distance(leaving_slots.begin(), leaving_end);
  leaving_slots.resize(static_cast<std::size_t>(leaving_particles));

  thrust::stable_sort(
    w.on_this(),
    leaving_slots.begin(),
    leaving_slots.end(),
    [=](const auto lhs, const auto rhs) {
      return slot_to_subregion_index(lhs) < slot_to_subregion_index(rhs);
    });

  // Gather their states and payloads into the buffers, then free their slots.
  const auto prev_buffer_size = buffer.size();
  buffer.resize(prev_buffer_size + leaving_slots.size());
  std::ignore = thrust::transform(
    w.on_this(),
    leaving_slots.begin(),
    leaving_slots.end(),
    rn::next(buffer.begin(), static_cast<std::ptrdiff_t>(prev_buffer_size)),
    [=](const auto n) {
      return runko::ParticleState<value_type> {
        .pos { pos_mds[n][0], pos_mds[n][1], pos_mds[n][2] },
        .vel { vel_mds[n][0], vel_mds[n][1], vel_mds[n][2] },
        .id { ids_mds[n][] },
        .wgt { weighted ? wgt_mds[n][] : w0 }
      };
    });

  // n_payload floats per leaving particle, in the order of the states
  const auto np = this->n_payload();
  const auto prev_payload_size = payload_buffer.size();
  payload_buffer.resize(prev_payload_size + leaving_slots.size() * np);
  const auto slots = thrust::raw_pointer_cast(leaving_slots.data());
  const auto dst   = rn::next(
    thrust::raw_pointer_cast(payload_buffer.data()),
    static_cast<std::ptrdiff_t>(prev_payload_size));
  for(auto c = 0uz; c < np; ++c) {
    thrust::for_each(
      w.on_this(),
      particle_ordinals_begin,
      rn::next(particle_ordinals_begin, leaving_particles),
      [=, src = this->payload_mds(c)](const auto k) {
        dst[k * np + c] = src[slots[k]][];
      });
  }

  thrust::for_each(
    w.on_this(),
    leaving_slots.begin(),
    leaving_slots.end(),
    [=](const auto n) { ids_mds[n][] = runko::dead_prtc_id; });

  using subregion_counter     = toolbox::VecD<runko::index_t, 27>;
  const auto slot_to_counter = [=](const runko::index_t n) {
    auto v                        = subregion_counter {};
    v[slot_to_subregion_index(n)] = 1;
    return v;
  };

  const auto subregion_counts = thrust::reduce(
    w.on_this(),
    thrust::transform_iterator(leaving_slots.begin(), slot_to_counter),
    thrust::transform_iterator(leaving_slots.end(), slot_to_counter),
    subregion_counter {});

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
  const auto kinetic_energies_end = rn::next(kinetic_energies_begin, this->ssize());

  return thrust::reduce(w.on_this(), kinetic_energies_begin, kinetic_energies_end);
}

ParticleContainer::ParticleContainer(const ParticleContainer& other)
{
  tyvi::copy_d2d(other.pos_, this->pos_);
  tyvi::copy_d2d(other.vel_, this->vel_);
  tyvi::copy_d2d(other.ids_, this->ids_);
  if(other.has_particle_weights()) { tyvi::copy_d2d(other.wgt_, this->wgt_); }
  this->payload_.resize(other.n_payload());
  for(auto c = 0uz; c < other.n_payload(); ++c) {
    tyvi::copy_d2d(other.payload_[c], this->payload_[c]);
  }

  this->charge_               = other.charge_;
  this->mass_                 = other.mass_;
  this->weight_               = other.weight_;
  this->reorder_count_ = other.reorder_count_;
}


ParticleContainer&
  ParticleContainer::operator=(const ParticleContainer& other)
{
  tyvi::copy_d2d(other.pos_, this->pos_);
  tyvi::copy_d2d(other.vel_, this->vel_);
  tyvi::copy_d2d(other.ids_, this->ids_);
  tyvi::copy_d2d(other.wgt_, this->wgt_);  // empty for an unweighted species
  this->payload_.resize(other.n_payload());
  for(auto c = 0uz; c < other.n_payload(); ++c) {
    tyvi::copy_d2d(other.payload_[c], this->payload_[c]);
  }

  this->charge_               = other.charge_;
  this->mass_                 = other.mass_;
  this->weight_               = other.weight_;
  this->reorder_count_ = other.reorder_count_;

  return *this;
}

}  // namespace pic
