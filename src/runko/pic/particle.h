// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/communication_common.h"
#include "runko/emf/yee_lattice.h"
#include "runko/mdgrid_common.h"
#include "runko/particles_common.h"
#include "runko/pic/reflector_wall.h"
#include "thrust/copy.h"
#include "thrust/count.h"
#include "thrust/device_vector.h"
#include "thrust/find.h"
#include "thrust/gather.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <any>
#include <array>
#include <concepts>
#include <cstddef>
#include <format>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <source_location>
#include <span>
#include <stdexcept>
#include <thrust/iterator/reverse_iterator.h>
#include <tuple>
#include <type_traits>
#include <utility>


namespace pic {

struct ParticleContainerArgs {
  std::size_t N;
  double charge, mass;
};

template<typename F, typename T>
concept score_function = std::regular_invocable<F, T, T, T> and
                         std::sortable<std::invoke_result_t<F, T, T, T>*>;

class [[nodiscard]] ParticleContainer {
public:
  /// The type in which pos and vel are stored in.
  using value_type = float;

  /// Specifies some span of particles in ParticleContainer.
  struct span {
    std::size_t begin { 0 }, end { 0 };
    [[nodiscard]] constexpr std::size_t size() const { return end - begin; }
  };

  /// Specifies a specifc span of particles in some ParticleContainer.
  ///
  /// Has to forward declared, because ParticleContainer
  /// is incomplete type at this point.
  struct specific_span;

private:
  /// Positions are in code units (i.e. x~ in x = x~ * Delta x).
  runko::device_vec_segments<value_type> pos_;
  /// Three-velocities are stored in physical/natural units.
  runko::device_vec_segments<value_type> vel_;
  /// Id for each particle.
  runko::device_scalar_segments<runko::prtc_id_type> ids_;

  double charge_;
  double mass_;

  /// Reusable buffer for particle trackers.
  thrust::device_vector<runko::index_t> tracker_cache_;
  /// Reusable buffer for generic data.
  thrust::device_vector<std::byte> generic_cache_;
  std::any sorting_score_cache_;

public:
  explicit ParticleContainer(ParticleContainerArgs);

  // Because tyvi::mdsegments is not copyable, we have to implement
  // copy ctor and copy assignment manually.

  ~ParticleContainer()                                       = default;
  explicit ParticleContainer(ParticleContainer&&) noexcept   = default;
  ParticleContainer& operator=(ParticleContainer&&) noexcept = default;

  explicit ParticleContainer(const ParticleContainer& other);
  ParticleContainer& operator=(const ParticleContainer& other);


  /// Construct the container by aggregating the given spans.
  template<std::ranges::forward_range R>
    requires std::convertible_to<
      std::ranges::range_reference_t<R>,
      ParticleContainer::specific_span>
  explicit ParticleContainer(const R&);

  /// Overwrite contents of this container with aggregated contents of the given spans.
  template<std::ranges::forward_range R>
    requires std::convertible_to<
      std::ranges::range_reference_t<R>,
      ParticleContainer::specific_span>
  void set(const R&);


  /// Appends the particles pointed by the given spans.
  template<std::ranges::forward_range R>
    requires std::convertible_to<
      std::ranges::range_reference_t<R>,
      std::span<const runko::ParticleState<value_type>>>
  inline void append(
    const R& spans,
    std::optional<std::array<value_type, 3>> wrap_mins = {},
    std::optional<std::array<value_type, 3>> wrap_maxs = {});

  /// Pre-allocate N dead-particle slots on an empty container.
  ///
  /// Sizes pos_/vel_/ids_ to N and tags every id with runko::dead_prtc_id.
  /// pos_ and vel_ are left uninitialized — dead particles are identified
  /// by their id alone and their position/velocity are never read. Must be
  /// called on a container where size() == 0. No temporary buffer is
  /// allocated, so peak memory during this call equals the final size.
  void prealloc_dead(std::size_t N);

  /// Returns the number of particles dead or alive.
  std::size_t size() const;

  /// Returns the number of particles dead or alive.
  std::ptrdiff_t ssize() const;

  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();
  std::vector<runko::prtc_id_type> get_ids();

  [[nodiscard]]
  auto pos_mds() const
  {
    return pos_.mds();
  }

  [[nodiscard]]
  auto vel_mds() const
  {
    return vel_.mds();
  }

  [[nodiscard]]
  auto ids_mds() const
  {
    return ids_.mds();
  }

  /// Wraps coordinates outside of given extents.
  void wrap_positions(std::array<value_type, 3> mins, std::array<value_type, 3> maxs);

  /// Splits container to 27 subcontainers along given divider lines.
  ///
  /// Data of the particles (x, y, z, vx, vy, vz, id)
  /// is appended to the given buffer in AOS manner.
  ///
  /// Returns map: (i, j, k) -> span
  /// where span represent the data of subregion (i, j, k) in the buffer.
  [[nodiscard]]

  std::map<runko::grid_neighbor<3>, std::pair<std::size_t, std::size_t>>
    divide_to_subregions(
      thrust::device_vector<runko::ParticleState<value_type>>& buffer,
      std::array<value_type, 2> x_dividers,
      std::array<value_type, 2> y_dividers,
      std::array<value_type, 2> z_dividers);

  /// Way to initialize particles. Known to be slow.
  template<std::ranges::forward_range R>
    requires std::
      convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState<double>>
    void add_particles(const R& new_particles);

  /// Sorts particles based on a function of particle position.
  ///
  /// Given function has to be device callable.
  void sort(pic::score_function<value_type> auto&& f);

  /// Push particles velocities and positions using boris scheme.
  inline void push_particles_boris(
    double cfl,
    runko::EB_interpolator<value_type> auto interpolator);

  /// Push particles velocities and positions using Higuera-Cary scheme.
  inline void push_particles_higuera_cary(
    double cfl,
    runko::EB_interpolator<value_type> auto interpolator);

  /// Push particles velocities and positions using Faraday-Cayley scheme.
  inline void push_particles_faraday(
    double cfl,
    runko::EB_interpolator<value_type> auto interpolator);

  emf::YeeLattice::CurrentContributions current_zigzag_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const double cfl) const;

  /// Adds generated current to given grid.
  ///
  /// `lattice_origo_coordinates` is in coordinates of the particles,
  /// and the spacing between cells in `Jout` is 1.
  /// Assumes that all particle contributions are inside the lattice.
  void current_zigzag_1st(
    runko::VecGrid<emf::YeeLattice::value_type>& Jout,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const double cfl) const;


  /// Reflect particles crossing a wall and deposit correction currents.
  ///
  /// Particles with global x-position behind walloc are reflected.
  /// Correction currents (both +q and -q zigzag deposits) are atomically
  /// added to correction_J so that the net current behind the wall is zero
  /// after the normal deposit_current() call.
  void reflect_at_wall(
    const reflector_wall& wall,
    runko::VecGrid<emf::YeeLattice::value_type>& correction_J,
    const std::array<value_type, 3> lattice_origo_coordinates,
    double cfl);

  /// Returns total kinetic energy of particless.
  ///
  /// Kinetic energy is given in units of mc^2,
  /// where m is the mass of the particle.
  double total_kinetic_energy() const;
};

struct ParticleContainer::specific_span {
  ParticleContainer::span span {};
  ParticleContainer const* container { nullptr };

  [[nodiscard]] constexpr std::size_t size() const { return span.size(); }
};


template<std::ranges::forward_range R>
  requires std::
    convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState<double>>
  inline void
  ParticleContainer::add_particles(const R& new_particles)
{
  using value_state_type = runko::ParticleState<ParticleContainer::value_type>;
  auto casted_new_particles =
    new_particles | std::views::transform([](const auto& state) {
      return value_state_type { .pos { static_cast<value_type>(state.pos[0]),
                                       static_cast<value_type>(state.pos[1]),
                                       static_cast<value_type>(state.pos[2]) },
                                .vel { static_cast<value_type>(state.vel[0]),
                                       static_cast<value_type>(state.vel[1]),
                                       static_cast<value_type>(state.vel[2]) },
                                .id = state.id };
    });

  const auto new_states_h = thrust::host_vector<value_state_type>(
    casted_new_particles.begin(),
    casted_new_particles.end());
  auto new_states = thrust::device_vector<value_state_type>(new_particles.size());
  thrust::copy(
    tyvi::mdgrid_work().on_this(),
    new_states_h.begin(),
    new_states_h.end(),
    new_states.begin());


  this->append(
    std::array { std::span<const value_state_type>(
      thrust::raw_pointer_cast(new_states.data()),
      new_states.size()) });
}

template<std::ranges::forward_range R>
  requires std::
    convertible_to<std::ranges::range_reference_t<R>, ParticleContainer::specific_span>
  inline ParticleContainer::ParticleContainer(const R& spans)
{
  this->set(spans);
}

template<std::ranges::forward_range R>
  requires std::
    convertible_to<std::ranges::range_reference_t<R>, ParticleContainer::specific_span>
  inline void
  ParticleContainer::set(const R& spans)
{
  namespace rn = std::ranges;
  namespace rv = std::views;

  if(rn::empty(spans)) {
    // Can not construct the container fron zero spans,
    // because there is no charge and mass specified.
    throw std::runtime_error { "Trying to set ParticleContainer from zero spans." };
  }

  /* Due to limitation of tyvi::when_all taking only statically known
     amount of tyvi::mdgrid_work, we have to wait for the copies in batches.
     (On a deeper level this is limitation of thrust.)

     However, this ctor takes dynamic amount of spans.
     To call tyvi::when_all on all of them, we have to do a ugly hack.
     In order to do so we can only support up to Nmax spans. */

  static constexpr auto Nmax = 27uz;  // atm nothing should call this with more.
  const auto Nspans          = rn::distance(spans);

  if(static_cast<std::size_t>(Nspans) > Nmax) {
    throw std::runtime_error {
      "ParticleContainer: trying to set from too many spans."
    };
  }

  {
    const auto args = rn::begin(spans)->container->args();
    this->charge_   = args.charge;
    this->mass_     = args.mass;

    const auto consistent_given_charge_n_mass =
      rn::all_of(spans, [this](const auto& span) {
        const auto args = span.container->args();
        return this->charge_ == args.charge and this->mass_ == args.mass;
      });
    if(not consistent_given_charge_n_mass) {
      throw std::runtime_error {
        "Trying to set ParticleContainer from spans to containers with "
        "inconsistent charges/masses."
      };
    }
  }

  const auto span_sizes =
    spans | rv::transform(&ParticleContainer::specific_span::size);

  // Where spans of particles begin and end in the resulting container:

  // { 0, span_sizes[0], ..., span_sizes[0] + ... + span_sizes[-2] }
  // Note that the last sum does not include span_sizes[-1].
  const auto begins = [&] {
    auto temp = std::vector<std::size_t>(Nspans);
    std::exclusive_scan(span_sizes.begin(), span_sizes.end(), temp.begin(), 0uz);
    return temp;
  }();

  // { span_sizes[0], ..., span_sizes[0] + ... + span_sizes[-1] }
  // Note that the last sum includes span_sizes[-1].
  const auto ends = [&] {
    auto temp = std::vector<std::size_t>(Nspans);
    std::inclusive_scan(span_sizes.begin(), span_sizes.end(), temp.begin());
    return temp;
  }();

  const auto Ntotal = ends.back();

  this->pos_.resize(Ntotal);
  this->vel_.resize(Ntotal);
  this->ids_.resize(Ntotal);

  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();
  const auto ids_mds = this->ids_.mds();

  const auto handle_span = [&](
                             tyvi::mdgrid_work& w,
                             const ParticleContainer::specific_span& span,
                             const std::array<std::size_t, 2> location_in_this) {
    const auto other_pos_mds = span.container->pos_.mds();
    const auto other_vel_mds = span.container->vel_.mds();
    const auto other_ids_mds = span.container->ids_.mds();

    const auto s                = std::array { span.span.begin, span.span.end };
    const auto other_pos_submds = std::submdspan(other_pos_mds, s);
    const auto other_vel_submds = std::submdspan(other_vel_mds, s);
    const auto other_ids_submds = std::submdspan(other_ids_mds, s);

    const auto pos_submds = std::submdspan(pos_mds, location_in_this);
    const auto vel_submds = std::submdspan(vel_mds, location_in_this);
    const auto ids_submds = std::submdspan(ids_mds, location_in_this);

    w.for_each_index(pos_submds, [=](const auto idx, const auto tidx) {
      pos_submds[idx][tidx] = other_pos_submds[idx][tidx];
      vel_submds[idx][tidx] = other_vel_submds[idx][tidx];
    });
    w.for_each_index(ids_submds, [=](const auto idx) {
      ids_submds[idx][] = other_ids_submds[idx][];
    });
  };

  auto works = std::vector<tyvi::mdgrid_work>(begins.size());
  for(auto i = 0uz; i < begins.size(); ++i) {
    handle_span(works[i], spans[i], { begins[i], ends[i] });
  }

  // Below is the ugly hack.
  // See beginning of the function for explanation why it is needed.

  auto hack = [&]<std::size_t... I>(std::index_sequence<I...>) {
    tyvi::when_all(works[I]...);
    // We assumed that there is at least one.
    works.front().wait();
  };

  switch(Nspans) {
    case 1: hack(std::make_index_sequence<1>()); break;
    case 2: hack(std::make_index_sequence<2>()); break;
    case 3: hack(std::make_index_sequence<3>()); break;
    case 4: hack(std::make_index_sequence<4>()); break;
    case 5: hack(std::make_index_sequence<5>()); break;
    case 6: hack(std::make_index_sequence<6>()); break;
    case 7: hack(std::make_index_sequence<7>()); break;
    case 8: hack(std::make_index_sequence<8>()); break;
    case 9: hack(std::make_index_sequence<9>()); break;
    case 10: hack(std::make_index_sequence<10>()); break;
    case 11: hack(std::make_index_sequence<11>()); break;
    case 12: hack(std::make_index_sequence<12>()); break;
    case 13: hack(std::make_index_sequence<13>()); break;
    case 14: hack(std::make_index_sequence<14>()); break;
    case 15: hack(std::make_index_sequence<15>()); break;
    case 16: hack(std::make_index_sequence<16>()); break;
    case 17: hack(std::make_index_sequence<17>()); break;
    case 18: hack(std::make_index_sequence<18>()); break;
    case 19: hack(std::make_index_sequence<19>()); break;
    case 20: hack(std::make_index_sequence<20>()); break;
    case 21: hack(std::make_index_sequence<21>()); break;
    case 22: hack(std::make_index_sequence<22>()); break;
    case 23: hack(std::make_index_sequence<23>()); break;
    case 24: hack(std::make_index_sequence<24>()); break;
    case 25: hack(std::make_index_sequence<25>()); break;
    case 26: hack(std::make_index_sequence<26>()); break;
    case 27: hack(std::make_index_sequence<27>()); break;
  };
}

template<std::ranges::forward_range R>
  requires std::convertible_to<
    std::ranges::range_reference_t<R>,
    std::span<const runko::ParticleState<ParticleContainer::value_type>>>
inline void
  ParticleContainer::append(
    const R& spans,
    const std::optional<std::array<value_type, 3>> wrap_mins,
    const std::optional<std::array<value_type, 3>> wrap_maxs)
{
  namespace rn = std::ranges;
  namespace rv = std::views;

  if(rn::empty(spans)) { return; }

  // 1. find last alive particle P
  // 2. resizes containers to fit appended particles (with tyvi::mdsegments this should
  // be cheap)
  // 3. copy appended particles to dead particles after P

  // 1.

  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->ssize());

  const auto rev_particle_ordinals_begin =
    thrust::reverse_iterator(particle_ordinals_end);
  const auto rev_particle_ordinals_end =
    thrust::reverse_iterator(particle_ordinals_begin);

  const auto w     = tyvi::mdgrid_work {};
  const auto Piter = thrust::find_if(
    w.on_this(),
    rev_particle_ordinals_begin,
    rev_particle_ordinals_end,
    [ids_mds = this->ids_.mds()](const auto n) {
      return ids_mds[n][] != runko::dead_prtc_id;
    });

  const auto dead_particles_at_end = rn::distance(rev_particle_ordinals_begin, Piter);
  const auto Pnum = static_cast<std::size_t>(this->ssize() - dead_particles_at_end);

  const auto span_sizes = spans | rv::transform(rn::size);

  // Where spans of particles begin P:
  // { 0, span_sizes[0], ..., span_sizes[0] + ... + span_sizes[-2] }
  // Note that the last sum does not include span_sizes[-1].
  const auto begins = [&] {
    auto temp = std::vector<std::size_t>(rn::size(spans));
    std::exclusive_scan(span_sizes.begin(), span_sizes.end(), temp.begin(), 0uz);
    return temp;
  }();

  // We can assume that there is at least one span,
  // because this was checked at the begin.
  const auto Ntotal = begins.back() + spans.back().size();

  // 2.
  this->pos_.resize(Pnum + Ntotal);
  this->vel_.resize(Pnum + Ntotal);
  this->ids_.resize(Pnum + Ntotal);

  // 3.
  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();
  const auto ids_mds = this->ids_.mds();

  for(auto n = 0uz; n < spans.size(); ++n) {
    const auto n_span   = spans[n];
    const auto n_size   = std::ranges::ssize(n_span);
    const auto n_begin  = thrust::counting_iterator<runko::index_t>(0uz);
    const auto n_end    = rn::next(n_begin, n_size);
    const auto n_offset = Pnum + begins[n];

    if(wrap_mins and wrap_maxs) {
      const auto ax = wrap_mins.value()[0];
      const auto ay = wrap_mins.value()[1];
      const auto az = wrap_mins.value()[2];
      const auto bx = wrap_maxs.value()[0];
      const auto by = wrap_maxs.value()[1];
      const auto bz = wrap_maxs.value()[2];
      const auto Lx = bx - ax;
      const auto Ly = by - ay;
      const auto Lz = bz - az;

      thrust::for_each(w.on_this(), n_begin, n_end, [=](const auto i) {
        // This does not work on lumi, so we use workaround.
        // const auto state = n_span[i];
        const auto state = n_span.data()[i];

        const auto j  = i + n_offset;
        pos_mds[j][0] = (state.pos[0] < 0 ? bx : ax) + std::fmod(state.pos[0], Lx);
        pos_mds[j][1] = (state.pos[1] < 0 ? by : ay) + std::fmod(state.pos[1], Ly);
        pos_mds[j][2] = (state.pos[2] < 0 ? bz : az) + std::fmod(state.pos[2], Lz);

        vel_mds[j][0] = state.vel[0];
        vel_mds[j][1] = state.vel[1];
        vel_mds[j][2] = state.vel[2];

        ids_mds[j][] = state.id;
      });

    } else {
      thrust::for_each(w.on_this(), n_begin, n_end, [=](const auto i) {
        // This does not work on lumi, so we use workaround.
        // const auto state = n_span[i];
        const auto state = n_span.data()[i];

        const auto j  = i + n_offset;
        pos_mds[j][0] = state.pos[0];
        pos_mds[j][1] = state.pos[1];
        pos_mds[j][2] = state.pos[2];

        vel_mds[j][0] = state.vel[0];
        vel_mds[j][1] = state.vel[1];
        vel_mds[j][2] = state.vel[2];

        ids_mds[j][] = state.id;
      });
    }

    w.wait();
  }
}

inline void
  ParticleContainer::sort(pic::score_function<value_type> auto&& f)
{
  /* This is overly verbose due to bug in gcc:
     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=12526

     Without it we could write:
     const auto pos0 = this->pos_.component_span<0>(); */

  using I                    = runko::index_t;
  static constexpr auto zero = std::array { I { 0 } };
  static constexpr auto one  = std::array { I { 1 } };
  static constexpr auto two  = std::array { I { 2 } };

  const auto pos0 = this->pos_.component_view<zero>();
  const auto pos1 = this->pos_.component_view<one>();
  const auto pos2 = this->pos_.component_view<two>();
  const auto vel0 = this->vel_.component_view<zero>();
  const auto vel1 = this->vel_.component_view<one>();
  const auto vel2 = this->vel_.component_view<two>();
  const auto ids  = this->ids_.component_view<>();

  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->ssize());
  this->tracker_cache_.assign(particle_ordinals_begin, particle_ordinals_end);

  using score_type =
    std::invoke_result_t<decltype(f), value_type, value_type, value_type>;

  const auto pos_mds = this->pos_.mds();
  const auto ids_mds = this->ids_.mds();

  const auto scores_begin = thrust::make_transform_iterator(
    particle_ordinals_begin,
    [=, f = std::forward<decltype(f)>(f)](const std::size_t i) {
      if(ids_mds[i][] == runko::dead_prtc_id) {
        return std::numeric_limits<score_type>::max();
      }
      return f(pos_mds[i][0], pos_mds[i][1], pos_mds[i][2]);
    });
  const auto scores_end = rn::next(scores_begin, this->ssize());

  using scores_cache_type = thrust::device_vector<score_type>;
  if(typeid(scores_cache_type) != this->sorting_score_cache_.type()) {
    this->sorting_score_cache_ = scores_cache_type {};
  }

  auto& scores_cache = std::any_cast<scores_cache_type&>(this->sorting_score_cache_);
  scores_cache.assign(scores_begin, scores_end);

  auto w = tyvi::mdgrid_work {};
  thrust::sort_by_key(
    w.on_this(),
    scores_cache.begin(),
    scores_cache.end(),
    this->tracker_cache_.begin());

  this->generic_cache_.resize(sizeof(value_type) * this->size());

  auto assert_alignment = []<typename T>(
                            const T* p,
                            const std::source_location loc =
                              std::source_location::current()) {
    if(reinterpret_cast<std::uintptr_t>(p) % alignof(T)) {
      const auto str = std::format("Incorrect alignment at line: {}", loc.line());
      throw std::runtime_error { str };
    }
  };

  const auto p = reinterpret_cast<value_type*>(
    thrust::raw_pointer_cast(this->generic_cache_.data()));
  assert_alignment(p);

  std::ignore = thrust::copy(w.on_this(), pos0.begin(), pos0.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    pos0.begin());
  std::ignore = thrust::copy(w.on_this(), pos1.begin(), pos1.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    pos1.begin());
  std::ignore = thrust::copy(w.on_this(), pos2.begin(), pos2.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    pos2.begin());
  std::ignore = thrust::copy(w.on_this(), vel0.begin(), vel0.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    vel0.begin());
  std::ignore = thrust::copy(w.on_this(), vel1.begin(), vel1.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    vel1.begin());
  std::ignore = thrust::copy(w.on_this(), vel2.begin(), vel2.end(), p);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    p,
    vel2.begin());

  this->generic_cache_.resize(sizeof(runko::prtc_id_type) * this->size());
  const auto pt = reinterpret_cast<runko::prtc_id_type*>(
    thrust::raw_pointer_cast(this->generic_cache_.data()));
  assert_alignment(pt);
  std::ignore = thrust::copy(w.on_this(), ids.begin(), ids.end(), pt);
  std::ignore = thrust::gather(
    w.on_this(),
    this->tracker_cache_.begin(),
    this->tracker_cache_.end(),
    pt,
    ids.begin());
  w.wait();
}
}  // namespace pic

#include "runko/pic/particle_boris.h"
#include "runko/pic/particle_faraday.h"
#include "runko/pic/particle_higuera_cary.h"
