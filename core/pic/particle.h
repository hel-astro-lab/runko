#pragma once

#include "core/emf/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "core/particles_common.h"
#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/sort.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <execution>
#include <functional>
#include <iterator>
#include <map>
#include <numeric>
#include <ranges>
#include <stdexcept>
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
  using E = std::dextents<std::size_t, 1>;

  runko::VecList<value_type> pos_;
  runko::VecList<value_type> vel_;

  double charge_;
  double mass_;

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

  /// Returns the number of particles.
  std::size_t size() const;
  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();

  [[nodiscard]]
  auto span_pos() &;

  [[nodiscard]]
  auto span_pos() const&;

  [[nodiscard]]
  auto span_vel() &;

  [[nodiscard]]
  auto span_vel() const&;

  /// Wraps coordinates outside of given extents.
  void wrap_positions(std::array<value_type, 3> mins, std::array<value_type, 3> maxs);

  /// Splits container to 27 subcontainers along given divider lines.
  ///
  /// Returns: [TODO | write this]
  [[nodiscard]]
  std::pair<std::map<std::array<int, 3>, ParticleContainer::span>, ParticleContainer>
    divide_to_subregions(
      std::array<value_type, 2> x_dividers,
      std::array<value_type, 2> y_dividers,
      std::array<value_type, 2> z_dividers);

  /// Add particles from other containers.
  template<std::convertible_to<const ParticleContainer&>... T>
  void add_particles(const T&... others);

  /// Way to initialize particles. Known to be slow.
  template<std::ranges::forward_range R>
    requires std::
      convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState>
    void add_particles(const R& new_particles);

  /// Sorts particles based on a function of particle position.
  ///
  /// Given function has to be device callable.
  void sort(pic::score_function<value_type> auto&& f);

  using InterpolatedEB_function =
    std::function<emf::YeeLattice::InterpolatedEB(const runko::VecList<value_type>&)>;

  /// Push particles velocities and positions using boris scheme.
  void push_particles_boris(double cfl, const InterpolatedEB_function&);

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
    std::exclusive_scan(
      std::execution::unseq,
      Nothers.begin(),
      Nothers.end(),
      temp.begin(),
      Nprev);
    return temp;
  }();

  // { Nprev + Nothers[0], ..., Nprev + Nothers[0] + ... + Nothers[-1] }
  // Note that the last sum includes Nothers[-1].
  const auto other_ends = [&] {
    auto temp = std::array<std::size_t, sizeof...(others)> {};
    std::inclusive_scan(
      std::execution::unseq,
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

  const auto new_pos_mds = new_pos.mds();
  const auto new_vel_mds = new_vel.mds();

  const auto prev_pos_mds = this->pos_.mds();
  const auto prev_vel_mds = this->vel_.mds();

  auto wA = tyvi::mdgrid_work {};
  wA.for_each_index(prev_pos_mds, [=](const auto idx, const auto tidx) {
    new_pos_mds[idx][tidx] = prev_pos_mds[idx][tidx];
    new_vel_mds[idx][tidx] = prev_vel_mds[idx][tidx];
  });

  const auto handle_other = [&](
                              tyvi::mdgrid_work& w,
                              const ParticleContainer& other,
                              const std::array<std::size_t, 2> where_other_goes) {
    const auto other_pos_mds = other.pos_.mds();
    const auto other_vel_mds = other.vel_.mds();

    const auto other_in_new_pos_mds = std::submdspan(new_pos_mds, where_other_goes);
    const auto other_in_new_vel_mds = std::submdspan(new_vel_mds, where_other_goes);

    w.for_each_index(other_in_new_pos_mds, [=](const auto idx, const auto tidx) {
      other_in_new_pos_mds[idx][tidx] = other_pos_mds[idx][tidx];
      other_in_new_vel_mds[idx][tidx] = other_vel_mds[idx][tidx];
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
}


template<std::ranges::forward_range R>
  requires std::convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState>
inline void
  ParticleContainer::add_particles(const R& new_particles)
{
  const auto N = static_cast<std::size_t>(std::ranges::distance(new_particles));

  auto args  = this->args();
  args.N     = N;
  auto added = ParticleContainer(args);

  const auto added_pos_smds = added.pos_.staging_mds();
  const auto added_vel_smds = added.vel_.staging_mds();

  for(const auto [i, p]: std::views::enumerate(new_particles)) {
    for(const auto j: std::views::iota(0uz, 3uz)) {
      added_pos_smds[i][j] = p.pos[j];
      added_vel_smds[i][j] = p.vel[j];
    }
  }

  auto w1 = tyvi::mdgrid_work {};
  auto w2 = tyvi::mdgrid_work {};
  w1.sync_from_staging(added.pos_);
  w2.sync_from_staging(added.vel_);

  tyvi::when_all(w1, w2);
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

  if(Nspans > Nmax) {
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
    std::exclusive_scan(
      std::execution::unseq,
      span_sizes.begin(),
      span_sizes.end(),
      temp.begin(),
      0uz);
    return temp;
  }();

  // { span_sizes[0], ..., span_sizes[0] + ... + span_sizes[-1] }
  // Note that the last sum includes span_sizes[-1].
  const auto ends = [&] {
    auto temp = std::vector<std::size_t>(Nspans);
    std::inclusive_scan(
      std::execution::unseq,
      span_sizes.begin(),
      span_sizes.end(),
      temp.begin());
    return temp;
  }();

  const auto Ntotal = ends.back();

  this->pos_.invalidating_resize(Ntotal);
  this->vel_.invalidating_resize(Ntotal);

  const auto pos_mds = this->pos_.mds();
  const auto vel_mds = this->vel_.mds();

  const auto handle_span = [&](
                             tyvi::mdgrid_work& w,
                             const ParticleContainer::specific_span& span,
                             const std::array<std::size_t, 2> location_in_this) {
    const auto other_pos_mds = span.container->pos_.mds();
    const auto other_vel_mds = span.container->vel_.mds();

    const auto s                = std::array { span.span.begin, span.span.end };
    const auto other_pos_submds = std::submdspan(other_pos_mds, s);
    const auto other_vel_submds = std::submdspan(other_vel_mds, s);

    const auto pos_submds = std::submdspan(pos_mds, location_in_this);
    const auto vel_submds = std::submdspan(vel_mds, location_in_this);

    w.for_each_index(pos_submds, [=](const auto idx, const auto tidx) {
      pos_submds[idx][tidx] = other_pos_submds[idx][tidx];
      vel_submds[idx][tidx] = other_vel_submds[idx][tidx];
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


inline void
  ParticleContainer::sort(pic::score_function<value_type> auto&& f)
{
  namespace rn                       = std::ranges;
  const auto particle_ordinals_begin = thrust::counting_iterator<std::size_t>(0uz);
  const auto particle_ordinals_end   = rn::next(particle_ordinals_begin, this->size());

  const auto pos_mds = this->pos_.mds();

  const auto scores_begin = thrust::make_transform_iterator(
    particle_ordinals_begin,
    [=, f = std::forward<decltype(f)>(f)](const std::size_t i) {
      return f(pos_mds[i][0], pos_mds[i][1], pos_mds[i][2]);
    });
  const auto scores_end = rn::next(scores_begin, this->size());

  auto trackers =
    thrust::device_vector<std::size_t>(particle_ordinals_begin, particle_ordinals_end);

  using score_t = std::invoke_result_t<decltype(f), value_type, value_type, value_type>;
  auto scores   = thrust::device_vector<score_t>(scores_begin, scores_end);

  thrust::sort_by_key(scores.begin(), scores.end(), trackers.begin());

  auto tmp_pos           = this->pos_;
  auto tmp_vel           = this->vel_;
  const auto tmp_pos_mds = tmp_pos.mds();
  const auto tmp_vel_mds = tmp_vel.mds();

  const auto vel_mds = this->vel_.mds();

  tyvi::mdgrid_work {}
    .for_each_index(
      tmp_pos,
      [=, p = trackers.begin()](const auto idx, const auto tidx) {
        const std::size_t i = p[idx[0]];
        pos_mds[idx][tidx]  = tmp_pos_mds[i][tidx];
        vel_mds[idx][tidx]  = tmp_vel_mds[i][tidx];
      })
    .wait();
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


}  // namespace pic
