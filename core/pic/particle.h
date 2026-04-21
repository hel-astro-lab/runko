#pragma once

#include "core/communication_common.h"
#include "core/emf/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "core/particles_common.h"
#include "core/pic/reflector_wall.h"
#include "thrust/copy.h"
#include "thrust/count.h"
#include "thrust/device_vector.h"
#include "thrust/find.h"
#include "thrust/gather.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"
#include "tools/signum.h"
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
  runko::VecList<value_type> pos_;
  /// Three-velocities are stored in physical/natural units.
  runko::VecList<value_type> vel_;
  /// Id for each particle.
  runko::ScalarList<runko::prtc_id_type> ids_;

  double charge_;
  double mass_;

  /// Reusable buffer for particle trackers.
  thrust::device_vector<runko::index_t> tracker_cache_;
  /// Reusable buffer for generic data.
  thrust::device_vector<std::byte> generic_cache_;
  std::any sorting_score_cache_;

public:
  explicit ParticleContainer(ParticleContainerArgs);

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
  inline void append(const R& spans);

  /// Returns the number of particles dead or alive.
  std::size_t size() const;
  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();
  std::vector<runko::prtc_id_type> get_ids();

  [[nodiscard]]
  auto span_pos() &;

  [[nodiscard]]
  auto span_pos() const&;

  [[nodiscard]]
  auto span_vel() &;

  [[nodiscard]]
  auto span_vel() const&;

  [[nodiscard]]
  auto span_id() &;

  [[nodiscard]]
  auto span_id() const&;

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

  /// Add particles from other containers.
  template<std::convertible_to<const ParticleContainer&>... T>
  void add_particles(const T&... others);

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
    const value_type cfl) const;


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
    value_type cfl);

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

template<std::convertible_to<const ParticleContainer&>... T>
inline void
  ParticleContainer::add_particles(const T&... others)
{
  if(((&static_cast<const ParticleContainer&>(others) == this) or ...)) {
    throw std::runtime_error {
      "Trying to add particles to ParticleContainer from itself."
    };
  }

  const auto Nprev = this->size();
  const auto Nothers =
    std::array { static_cast<const ParticleContainer&>(others).size()... };

  // { Nprev, Nprev + Nothers[0], ..., Nprev + Nothers[0] + ... + Nothers[-2] }
  // Note that the last sum does not include Nothers[-1].
  const auto other_begins = [&] {
    auto temp = std::array<std::size_t, sizeof...(others)> {};
    std::exclusive_scan(Nothers.begin(), Nothers.end(), temp.begin(), Nprev);
    return temp;
  }();

  // { Nprev + Nothers[0], ..., Nprev + Nothers[0] + ... + Nothers[-1] }
  // Note that the last sum includes Nothers[-1].
  const auto other_ends = [&] {
    auto temp = std::array<std::size_t, sizeof...(others)> {};
    std::inclusive_scan(
      Nothers.begin(),
      Nothers.end(),
      temp.begin(),
      std::plus<> {},
      Nprev);
    return temp;
  }();

  const auto Ntotal = other_ends.back();

  auto new_pos = runko::VecList<value_type>(Ntotal);
  auto new_vel = runko::VecList<value_type>(Ntotal);
  auto new_ids = runko::ScalarList<runko::prtc_id_type>(Ntotal);

  const auto new_pos_mds = new_pos.mds();
  const auto new_vel_mds = new_vel.mds();
  const auto new_ids_mds = new_ids.mds();

  const auto prev_pos_mds = this->pos_.mds();
  const auto prev_vel_mds = this->vel_.mds();
  const auto prev_ids_mds = this->ids_.mds();

  auto wA = tyvi::mdgrid_work {};
  wA.for_each_index(prev_pos_mds, [=](const auto idx, const auto tidx) {
    new_pos_mds[idx][tidx] = prev_pos_mds[idx][tidx];
    new_vel_mds[idx][tidx] = prev_vel_mds[idx][tidx];
  });
  wA.for_each_index(prev_ids_mds, [=](const auto idx) {
    new_ids_mds[idx][] = prev_ids_mds[idx][];
  });

  const auto handle_other = [&](
                              tyvi::mdgrid_work& w,
                              const ParticleContainer& other,
                              const std::array<std::size_t, 2> where_other_goes) {
    const auto other_pos_mds = other.pos_.mds();
    const auto other_vel_mds = other.vel_.mds();
    const auto other_ids_mds = other.ids_.mds();

    const auto other_in_new_pos_mds = std::submdspan(new_pos_mds, where_other_goes);
    const auto other_in_new_vel_mds = std::submdspan(new_vel_mds, where_other_goes);
    const auto other_in_new_ids_mds = std::submdspan(new_ids_mds, where_other_goes);

    w.for_each_index(other_in_new_pos_mds, [=](const auto idx, const auto tidx) {
      other_in_new_pos_mds[idx][tidx] = other_pos_mds[idx][tidx];
      other_in_new_vel_mds[idx][tidx] = other_vel_mds[idx][tidx];
    });

    w.for_each_index(other_in_new_ids_mds, [=](const auto idx) {
      other_in_new_ids_mds[idx][] = other_ids_mds[idx][];
    });
  };

  // Ugly syntax until C++26.
  [&]<std::size_t... I>(std::index_sequence<I...>) {
    auto other_works = std::array { (std::ignore = I, tyvi::mdgrid_work {})... };
    (handle_other(other_works[I], others, { other_begins[I], other_ends[I] }), ...);
    tyvi::when_all(wA, other_works[I]...);
    wA.wait();
  }(std::make_index_sequence<sizeof...(others)>());

  pos_ = std::move(new_pos);
  vel_ = std::move(new_vel);
  ids_ = std::move(new_ids);
}


template<std::ranges::forward_range R>
  requires std::
    convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState<double>>
  inline void
  ParticleContainer::add_particles(const R& new_particles)
{
  const auto N = static_cast<std::size_t>(std::ranges::distance(new_particles));

  auto args  = this->args();
  args.N     = N;
  auto added = ParticleContainer(args);

  const auto added_pos_smds = added.pos_.staging_mds();
  const auto added_vel_smds = added.vel_.staging_mds();
  const auto added_ids_smds = added.ids_.staging_mds();

  {
    std::size_t i = 0;
    for(const auto& p: new_particles) {
      for(const auto j: std::views::iota(0uz, 3uz)) {
        added_pos_smds[i][j] = p.pos[j];
        added_vel_smds[i][j] = p.vel[j];
      }

      added_ids_smds[i][] = p.id;
      ++i;
    }
  }

  auto w1 = tyvi::mdgrid_work {};
  auto w2 = tyvi::mdgrid_work {};
  auto w3 = tyvi::mdgrid_work {};
  w1.sync_from_staging(added.pos_);
  w2.sync_from_staging(added.vel_);
  w3.sync_from_staging(added.ids_);

  tyvi::when_all(w1, w2, w3);
  w1.wait();

  this->add_particles(added);
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

  this->pos_.invalidating_resize(Ntotal);
  this->vel_.invalidating_resize(Ntotal);
  this->ids_.invalidating_resize(Ntotal);

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
  ParticleContainer::append(const R& spans)
{
  namespace rn = std::ranges;
  namespace rv = std::views;

  if(rn::empty(spans)) { return; }

  // 1. find last alive particle P
  // 2. check if the appended particles fit to dead particles after P
  // 3. if no, resize such that they fit
  // 4. copy appended particles to dead particles after P

  // 1.

  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->size());

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
  const auto Pnum                  = this->size() - dead_particles_at_end;

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
  if(Ntotal > static_cast<std::size_t>(dead_particles_at_end)) {
    // 3.

    // This could be more efficient, as now we copy all the staging buffers as well.
    const auto tmp = *this;

    const auto tmp_pos_mds = tmp.pos_.mds();
    const auto tmp_vel_mds = tmp.vel_.mds();
    const auto tmp_ids_mds = tmp.ids_.mds();

    this->pos_.invalidating_resize(Pnum + Ntotal);
    this->vel_.invalidating_resize(Pnum + Ntotal);
    this->ids_.invalidating_resize(Pnum + Ntotal);

    const auto pos_mds = this->pos_.mds();
    const auto vel_mds = this->vel_.mds();
    const auto ids_mds = this->ids_.mds();

    w.for_each_index(tmp_pos_mds, [=](const auto idx, const auto tidx) {
      pos_mds[idx][tidx] = tmp_pos_mds[idx][tidx];
      vel_mds[idx][tidx] = tmp_vel_mds[idx][tidx];
    });
    w.for_each_index(tmp_ids_mds, [=](const auto idx) {
      ids_mds[idx][] = tmp_ids_mds[idx][];
    });
    w.wait();
  }
  // 4.

  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();
  const auto ids_mds = this->ids_.mds();

  for(auto n = 0uz; n < spans.size(); ++n) {
    const auto n_span   = spans[n];
    const auto n_size   = n_span.size();
    const auto n_begin  = thrust::counting_iterator<runko::index_t>(0uz);
    const auto n_end    = rn::next(n_begin, n_size);
    const auto n_offset = Pnum + begins[n];

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
    w.wait();
  }
}

inline void
  ParticleContainer::sort(pic::score_function<value_type> auto&& f)
{
  const auto pos0 = this->pos_.device_component_span<0>();
  const auto pos1 = this->pos_.device_component_span<1>();
  const auto pos2 = this->pos_.device_component_span<2>();
  const auto vel0 = this->vel_.device_component_span<0>();
  const auto vel1 = this->vel_.device_component_span<1>();
  const auto vel2 = this->vel_.device_component_span<2>();
  const auto ids  = this->ids_.device_component_span<>();

  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<runko::index_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->size());
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
  const auto scores_end = rn::next(scores_begin, this->size());

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

inline auto
  ParticleContainer::span_pos() &
{
  return this->pos_.span();
}

inline auto
  ParticleContainer::span_pos() const&
{
  return this->pos_.span();
}

inline auto
  ParticleContainer::span_vel() &
{
  return this->vel_.span();
}

inline auto
  ParticleContainer::span_vel() const&
{
  return this->vel_.span();
}

inline auto
  ParticleContainer::span_id() &
{
  return this->ids_.span();
}

inline auto
  ParticleContainer::span_id() const&
{
  return this->ids_.span();
}

}  // namespace pic

#include "core/pic/particle_boris.h"
#include "core/pic/particle_faraday.h"
#include "core/pic/particle_higuera_cary.h"
