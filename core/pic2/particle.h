#pragma once

#include "core/emf2/yee_lattice.h"
#include "core/mdgrid_common.h"
#include "core/particles_common.h"
#include "thrust/device_vector.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <execution>
#include <functional>
#include <numeric>
#include <ranges>
#include <type_traits>

namespace pic2 {

struct ParticleContainerArgs {
  std::size_t N;
  double charge, mass;
};

class [[nodiscard]] ParticleContainer {
public:
  /// The type in which pos and vel are stored in.
  using value_type = float;

private:
  using E = std::dextents<std::size_t, 1>;

  runko::VecList<value_type> pos_;
  runko::VecList<value_type> vel_;

  double charge_;
  double mass_;

public:
  explicit ParticleContainer(ParticleContainerArgs);

  /// Returns the number of particles.
  std::size_t size() const;
  ParticleContainerArgs args() const;

  std::array<std::vector<value_type>, 3> get_positions();
  std::array<std::vector<value_type>, 3> get_velocities();

  [[nodiscard]]
  auto span_pos(this auto& self)
  {
    return self.pos_.span();
  }

  [[nodiscard]]
  auto span_vel(this auto& self)
  {
    return self.vel_.span();
  }

  /// Splits container to 27 subcontainers along given divider lines.
  ///
  /// Container pointed by this is replaced by the middle container
  /// and the rest are returned with corresponding direction.
  ///
  /// Applies periodic boundary condition based on global_{mins,maxs}.
  std::vector<std::pair<std::array<int, 3>, ParticleContainer>> split_to_subregions(
    std::array<value_type, 2> x_dividers,
    std::array<value_type, 2> y_dividers,
    std::array<value_type, 2> z_dividers,
    std::array<value_type, 3> global_mins,
    std::array<value_type, 3> global_maxs);

  /// Add particles from other containers.
  template<std::convertible_to<const ParticleContainer&>... T>
  void add_particles(const T&... others);

  /// Way to initialize particles. Known to be slow.
  template<std::ranges::forward_range R>
    requires std::
      convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState>
    void add_particles(const R& new_particles);

  using InterpolatedEB_function =
    std::function<emf2::YeeLattice::InterpolatedEB(const runko::VecList<value_type>&)>;

  /// Push particles velocities and positions using boris scheme.
  void push_particles_boris(double cfl, const InterpolatedEB_function&);

  emf2::YeeLattice::CurrentContributions current_zigzag_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const;

  /// Adds generated current to given grid.
  ///
  /// `lattice_orgio_coordinates` is in coordinates of the particles,
  /// and the spacing between cells in `Jout` is 1.
  /// Assumes that all particle contributions are inside the lattice.
  void current_zigzag_1st(
    runko::VecGrid<emf2::YeeLattice::value_type>& Jout,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const;
};

/// Way to initialize particles. Known to be slow.
template<std::ranges::forward_range R>
  requires std::convertible_to<std::ranges::range_reference_t<R>, runko::ParticleState>
void
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

  auto w1 = tyvi::mdgrid_work {}.sync_from_staging(added.pos_);
  auto w2 = tyvi::mdgrid_work {}.sync_from_staging(added.vel_);
  tyvi::when_all(w1, w2).wait();

  this->add_particles(added);
}

/* It would potentially be convinient for this to take
   a range of ParticleContainers but due to limitation
   of tyvi::when_all not handling dynamic number of work items,
   the number of ParticleContainers has to be known statically. */
template<std::convertible_to<const ParticleContainer&>... T>
void
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

  auto wA = tyvi::mdgrid_work {}.for_each_index(
    prev_pos_mds,
    [=](const auto idx, const auto tidx) {
      new_pos_mds[idx][tidx] = prev_pos_mds[idx][tidx];
      new_vel_mds[idx][tidx] = prev_vel_mds[idx][tidx];
    });

  const auto handle_other = [&](
                              const ParticleContainer& other,
                              const std::array<std::size_t, 2> where_other_goes) {
    const auto other_pos_mds = other.pos_.mds();
    const auto other_vel_mds = other.vel_.mds();

    const auto other_in_new_pos_mds = std::submdspan(new_pos_mds, where_other_goes);
    const auto other_in_new_vel_mds = std::submdspan(new_vel_mds, where_other_goes);

    return tyvi::mdgrid_work {}.for_each_index(
      other_in_new_pos_mds,
      [=](const auto idx, const auto tidx) {
        other_in_new_pos_mds[idx][tidx] = other_pos_mds[idx][tidx];
        other_in_new_vel_mds[idx][tidx] = other_vel_mds[idx][tidx];
      });
  };

  // Ugly syntax until C++26.
  [&]<std::size_t... I>(std::index_sequence<I...>) {
    auto other_works =
      std::array { handle_other(others, { other_begins[I], other_ends[I] })... };
    tyvi::when_all(wA, other_works[I]...).wait();
  }(std::make_index_sequence<sizeof...(others)>());

  pos_ = std::move(new_pos);
  vel_ = std::move(new_vel);
}


}  // namespace pic2
