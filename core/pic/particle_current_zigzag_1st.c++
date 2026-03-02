#include "core/emf/yee_lattice.h"
#include "core/pic/particle.h"
#include "tools/vector.h"
#include "tyvi/algo.h"
#include "tyvi/containers.h"
#include "tyvi/iterators.h"
#include "tyvi/mdgrid.h"
#include "tyvi/sstd.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <utility>


namespace pic {

emf::YeeLattice::CurrentContributions
  ParticleContainer::current_zigzag_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const double cfl_) const
{
  using Vec3uz = toolbox::Vec3<uint32_t>;
  using Vec3v  = toolbox::Vec3<value_type>;

  //--------------------------------------------------
  // handy tool to process deposit location
  struct deposit_loc {
    Vec3uz idx;

    // Can not be defaulted due to amdgcn linker error for missing memcmp??
    // constexpr bool operator==(const deposit_loc&) const = default;
    constexpr bool operator==(const deposit_loc& rhs) const
    {
      return this->idx == rhs.idx;
    };

    constexpr bool operator<(const deposit_loc& rhs) const
    {
      const auto [i, j, k] = this->idx.data;
      const auto [a, b, c] = rhs.idx.data;

      return (i < a) or (i == a and j < b) or (i == a and j == b and k < c);
    }
  };

  //--------------------------------------------------
  static constexpr uint32_t max_num_of_nodes_effected = 14u;
  const auto N = this->size() * max_num_of_nodes_effected;

  auto deposit_locations     = tyvi::device_vector<deposit_loc>(N);
  const auto deposit_loc_ptr = deposit_locations.begin();

  //--------------------------------------------------
  // get current and particle array; define physical multipliers
  auto currents          = tyvi::device_vector<Vec3v>(N);
  const auto current_ptr = currents.begin();

  const auto pos_mds      = pos_.mds();
  const auto vel_mds      = vel_.mds();
  const value_type charge = charge_;
  const value_type cfl    = cfl_;

  //--------------------------------------------------
  // main loop over particles
  auto w = tyvi::mdgrid_work {};
  w.for_each_index(pos_mds, [=](const auto idx) {
    const auto u = Vec3v(vel_mds[idx]).template as<value_type>();
    const value_type invgam =
      value_type { 1 } / tyvi::sstd::sqrt(value_type { 1 } + toolbox::dot(u, u));

    const auto x2 = Vec3v(pos_mds[idx]).template as<value_type>() -
                    Vec3v(lattice_origo_coordinates).template as<value_type>();

    const auto x1 = x2 - cfl * invgam * u;

    // Float floor for relay and weight computation (pure float — no 64-bit integers)
    const auto fi1 = Vec3v(
      tyvi::sstd::floor(x1(0)), tyvi::sstd::floor(x1(1)), tyvi::sstd::floor(x1(2)));
    const auto fi2 = Vec3v(
      tyvi::sstd::floor(x2(0)), tyvi::sstd::floor(x2(1)), tyvi::sstd::floor(x2(2)));

    const auto relay = [&](const std::size_t j) -> value_type {
      const auto a  = tyvi::sstd::min(fi1(j), fi2(j)) + value_type { 1 };
      const auto b1 = tyvi::sstd::max(fi1(j), fi2(j));
      const auto b2 = value_type { 0.5 } * (x1(j) + x2(j));
      const auto b  = tyvi::sstd::max(b1, b2);
      return tyvi::sstd::min(a, b);
    };

    const auto x_relay = Vec3v(relay(0), relay(1), relay(2));

    const auto F1 = charge * (x_relay - x1);
    const auto F2 = charge * (x2 - x_relay);

    const auto W1 = value_type { 0.5 } * (x1 + x_relay) - fi1;
    const auto W2 = value_type { 0.5 } * (x2 + x_relay) - fi2;

    const auto [Fx1, Fy1, Fz1] = F1.data;
    const auto [Fx2, Fy2, Fz2] = F2.data;
    const auto [Wx1, Wy1, Wz1] = W1.data;
    const auto [Wx2, Wy2, Wz2] = W2.data;

    // Integer indices only for store_current (uint32_t = 32-bit = SIMD width=4)
    const auto i1 = fi1.template as<uint32_t>();
    const auto i2 = fi2.template as<uint32_t>();

    auto offset = static_cast<uint32_t>(idx[0]) * max_num_of_nodes_effected;

    const auto store_current = [=, &offset](const Vec3uz index, const Vec3v current) {
      deposit_loc_ptr[offset] = deposit_loc { index };
      current_ptr[offset++]   = current;
    };

    store_current(
      i1,
      Vec3v(
        Fx1 * (1.0f - Wy1) * (1.0f - Wz1),
        Fy1 * (1.0f - Wx1) * (1.0f - Wz1),
        Fz1 * (1.0f - Wx1) * (1.0f - Wy1)));

    store_current(
      i2,
      Vec3v(
        Fx2 * (1.0f - Wy2) * (1.0f - Wz2),
        Fy2 * (1.0f - Wx2) * (1.0f - Wz2),
        Fz2 * (1.0f - Wx2) * (1.0f - Wy2)));

    store_current(
      i1 + Vec3uz(1, 0, 0),
      Vec3v(0, Fy1 * Wx1 * (1.0f - Wz1), Fz1 * Wx1 * (1.0f - Wy1)));

    store_current(
      i2 + Vec3uz(1, 0, 0),
      Vec3v(0, Fy2 * Wx2 * (1.0f - Wz2), Fz2 * Wx2 * (1.0f - Wy2)));

    store_current(
      i1 + Vec3uz(0, 1, 0),
      Vec3v(Fx1 * Wy1 * (1.0f - Wz1), 0, Fz1 * (1.0f - Wx1) * Wy1));

    store_current(
      i2 + Vec3uz(0, 1, 0),
      Vec3v(Fx2 * Wy2 * (1.0f - Wz2), 0, Fz2 * (1.0f - Wx2) * Wy2));

    store_current(
      i1 + Vec3uz(0, 0, 1),
      Vec3v(Fx1 * (1.0f - Wy1) * Wz1, Fy1 * (1.0f - Wx1) * Wz1, 0));

    store_current(
      i2 + Vec3uz(0, 0, 1),
      Vec3v(Fx2 * (1.0f - Wy2) * Wz2, Fy2 * (1.0f - Wx2) * Wz2, 0));

    store_current(i1 + Vec3uz(0, 1, 1), Vec3v(Fx1 * Wy1 * Wz1, 0, 0));

    store_current(i2 + Vec3uz(0, 1, 1), Vec3v(Fx2 * Wy2 * Wz2, 0, 0));

    store_current(i1 + Vec3uz(1, 0, 1), Vec3v(0, Fy1 * Wx1 * Wz1, 0));

    store_current(i2 + Vec3uz(1, 0, 1), Vec3v(0, Fy2 * Wx2 * Wz2, 0));

    store_current(i1 + Vec3uz(1, 1, 0), Vec3v(0, 0, Fz1 * Wx1 * Wy1));

    store_current(i2 + Vec3uz(1, 1, 0), Vec3v(0, 0, Fz2 * Wx2 * Wy2));
  });

  auto unique_deposit_locations = tyvi::device_vector<std::array<std::size_t, 3>>(N);
  auto reduced_currents         = tyvi::device_vector<std::array<value_type, 3>>(N);

  const auto uniq_dep_locs_ptr = tyvi::make_transform_output_iterator(
    unique_deposit_locations.begin(),
    [](const deposit_loc& loc) -> std::array<std::size_t, 3> {
      return { loc.idx(0), loc.idx(1), loc.idx(2) };
    });

  const auto reduced_currents_ptr = tyvi::make_transform_output_iterator(
    reduced_currents.begin(),
    [](const Vec3v& J) { return J.data; });


  w.wait();

  tyvi::algo::sort_by_key(
    deposit_locations.begin(),
    deposit_locations.end(),
    currents.begin());

  const auto [key_end, _] = tyvi::algo::reduce_by_key(
    deposit_locations.begin(),
    deposit_locations.end(),
    currents.begin(),
    uniq_dep_locs_ptr,
    reduced_currents_ptr,
    std::equal_to<deposit_loc> {},
    std::plus<toolbox::Vec3<value_type>> {});

  const auto num_of_uniq_locs = key_end - uniq_dep_locs_ptr;

  unique_deposit_locations.resize(num_of_uniq_locs);
  reduced_currents.resize(num_of_uniq_locs);

  return emf::YeeLattice::CurrentContributions { .locations =
                                                   std::move(unique_deposit_locations),
                                                 .currents =
                                                   std::move(reduced_currents) };
}

void
  ParticleContainer::current_zigzag_1st(
    runko::VecGrid<emf::YeeLattice::value_type>& Jout,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const
{
  using Vec3uz = toolbox::Vec3<uint32_t>;
  using Vec3v  = toolbox::Vec3<value_type>;

  const auto Jmds         = Jout.mds();
  const auto pos_mds      = pos_.mds();
  const auto vel_mds      = vel_.mds();
  const value_type charge = charge_;

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        const auto u = Vec3v(vel_mds[idx]).template as<value_type>();
        const auto invgam =
          value_type { 1 } / tyvi::sstd::sqrt(value_type { 1 } + toolbox::dot(u, u));

        const auto x2 = Vec3v(pos_mds[idx]).template as<value_type>() -
                        Vec3v(lattice_origo_coordinates).template as<value_type>();

        const auto x1 = x2 - cfl * invgam * u;

        // Float floor for relay and weight computation (pure float — no 64-bit integers)
        const auto fi1 = Vec3v(
          tyvi::sstd::floor(x1(0)), tyvi::sstd::floor(x1(1)), tyvi::sstd::floor(x1(2)));
        const auto fi2 = Vec3v(
          tyvi::sstd::floor(x2(0)), tyvi::sstd::floor(x2(1)), tyvi::sstd::floor(x2(2)));

        const auto relay = [&](const std::size_t j) -> value_type {
          const auto a  = tyvi::sstd::min(fi1(j), fi2(j)) + value_type { 1 };
          const auto b1 = tyvi::sstd::max(fi1(j), fi2(j));
          const auto b2 = value_type { 0.5 } * (x1(j) + x2(j));
          const auto b  = tyvi::sstd::max(b1, b2);
          return tyvi::sstd::min(a, b);
        };

        const auto x_relay = Vec3v(relay(0), relay(1), relay(2));

        const auto F1 = charge * (x_relay - x1);
        const auto F2 = charge * (x2 - x_relay);

        const auto W1 = value_type { 0.5 } * (x1 + x_relay) - fi1;
        const auto W2 = value_type { 0.5 } * (x2 + x_relay) - fi2;

        const auto [Fx1, Fy1, Fz1] = F1.data;
        const auto [Fx2, Fy2, Fz2] = F2.data;
        const auto [Wx1, Wy1, Wz1] = W1.data;
        const auto [Wx2, Wy2, Wz2] = W2.data;

        // Integer indices only for store_current (uint32_t = 32-bit = SIMD width=4)
        const auto i1 = fi1.template as<uint32_t>();
        const auto i2 = fi2.template as<uint32_t>();

        const auto store_current = [&](const Vec3uz& index, const Vec3v& current) {
          const auto si  = index.template as<std::size_t>();
          auto* const Jx = &tyvi::sstd::raw_ref(Jmds[si.data][0]);
          auto* const Jy = &tyvi::sstd::raw_ref(Jmds[si.data][1]);
          auto* const Jz = &tyvi::sstd::raw_ref(Jmds[si.data][2]);

          tyvi::sstd::atomic_add(Jx, current[0]);
          tyvi::sstd::atomic_add(Jy, current[1]);
          tyvi::sstd::atomic_add(Jz, current[2]);
        };

        store_current(
          i1,
          Vec3v(
            Fx1 * (1.0f - Wy1) * (1.0f - Wz1),
            Fy1 * (1.0f - Wx1) * (1.0f - Wz1),
            Fz1 * (1.0f - Wx1) * (1.0f - Wy1)));

        store_current(
          i2,
          Vec3v(
            Fx2 * (1.0f - Wy2) * (1.0f - Wz2),
            Fy2 * (1.0f - Wx2) * (1.0f - Wz2),
            Fz2 * (1.0f - Wx2) * (1.0f - Wy2)));

        store_current(
          i1 + Vec3uz(1, 0, 0),
          Vec3v(0, Fy1 * Wx1 * (1.0f - Wz1), Fz1 * Wx1 * (1.0f - Wy1)));

        store_current(
          i2 + Vec3uz(1, 0, 0),
          Vec3v(0, Fy2 * Wx2 * (1.0f - Wz2), Fz2 * Wx2 * (1.0f - Wy2)));

        store_current(
          i1 + Vec3uz(0, 1, 0),
          Vec3v(Fx1 * Wy1 * (1.0f - Wz1), 0, Fz1 * (1.0f - Wx1) * Wy1));

        store_current(
          i2 + Vec3uz(0, 1, 0),
          Vec3v(Fx2 * Wy2 * (1.0f - Wz2), 0, Fz2 * (1.0f - Wx2) * Wy2));

        store_current(
          i1 + Vec3uz(0, 0, 1),
          Vec3v(Fx1 * (1.0f - Wy1) * Wz1, Fy1 * (1.0f - Wx1) * Wz1, 0));

        store_current(
          i2 + Vec3uz(0, 0, 1),
          Vec3v(Fx2 * (1.0f - Wy2) * Wz2, Fy2 * (1.0f - Wx2) * Wz2, 0));

        store_current(i1 + Vec3uz(0, 1, 1), Vec3v(Fx1 * Wy1 * Wz1, 0, 0));

        store_current(i2 + Vec3uz(0, 1, 1), Vec3v(Fx2 * Wy2 * Wz2, 0, 0));

        store_current(i1 + Vec3uz(1, 0, 1), Vec3v(0, Fy1 * Wx1 * Wz1, 0));

        store_current(i2 + Vec3uz(1, 0, 1), Vec3v(0, Fy2 * Wx2 * Wz2, 0));

        store_current(i1 + Vec3uz(1, 1, 0), Vec3v(0, 0, Fz1 * Wx1 * Wy1));

        store_current(i2 + Vec3uz(1, 1, 0), Vec3v(0, 0, Fz2 * Wx2 * Wy2));
      })
    .wait();
}

}  // namespace pic
