#include "core/emf2/yee_lattice.h"
#include "core/pic2/particle.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/iterator/transform_output_iterator.h"
#include "thrust/memory.h"
#include "thrust/reduce.h"
#include "thrust/sort.h"
#include "tools/vector.h"
#include "tyvi/mdgrid.h"

#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <print>
#include <stdexcept>
#include <utility>


namespace pic2 {

emf2::YeeLattice::CurrentContributions
  ParticleContainer::current_zigzag_1st(
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const
{
  using Vec3uz         = toolbox::Vec3<std::size_t>;
  using yee_value_type = emf2::YeeLattice::value_type;
  using Vec3yee        = toolbox::Vec3<yee_value_type>;
  using Vec3val        = toolbox::Vec3<value_type>;
  using Vec3d          = toolbox::Vec3<double>;

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

  static constexpr auto max_num_of_nodes_effected = 14uz;
  const auto N = this->size() * max_num_of_nodes_effected;

  auto deposit_locations     = thrust::device_vector<deposit_loc>(N);
  const auto deposit_loc_ptr = deposit_locations.begin();

  auto currents          = thrust::device_vector<Vec3yee>(N);
  const auto current_ptr = currents.begin();

  const auto pos_mds = pos_.mds();
  const auto vel_mds = vel_.mds();
  const auto charge  = charge_;

  auto w1 = tyvi::mdgrid_work {}.for_each_index(pos_mds, [=](const auto idx) {
    // NOTE: performing velocity calculations via doubles to retain accuracy

    const auto u      = Vec3val(vel_mds[idx]).template as<double>();
    const auto invgam = 1.0 / std::sqrt(1.0 + toolbox::dot(u, u));

    const auto x2 = Vec3val(pos_mds[idx]).template as<double>() -
                    Vec3val(lattice_origo_coordinates).template as<double>();

    const auto x1 = x2 - cfl * invgam * u;

    // These are equivalent to floor assuming that x1 and x2 do not have negative
    // values
    const auto i1 = x1.template as<std::size_t>();
    const auto i2 = x2.template as<std::size_t>();

    const auto relay = [&](const std::size_t j) -> double {
      const auto a  = static_cast<double>(min(i1(j), i2(j)) + 1);
      const auto b1 = static_cast<double>(max(i1(j), i2(j)));
      const auto b2 = 0.5 * (x1(j) + x2(j));
      const auto b  = max(b1, b2);
      return min(a, b);
    };

    const auto x_relay = Vec3d(relay(0), relay(1), relay(2));

    const auto F1 = charge * (x_relay - x1);
    const auto F2 = charge * (x2 - x_relay);

    const auto W1 = 0.5 * (x1 + x_relay) - i1.template as<double>();
    const auto W2 = 0.5 * (x2 + x_relay) - i2.template as<double>();

    const auto [Fx1, Fy1, Fz1] = F1.data;
    const auto [Fx2, Fy2, Fz2] = F2.data;
    const auto [Wx1, Wy1, Wz1] = W1.data;
    const auto [Wx2, Wy2, Wz2] = W2.data;

    auto offset = idx[0] * max_num_of_nodes_effected;

    const auto store_current = [=, &offset](const Vec3uz index, const Vec3yee current) {
      deposit_loc_ptr[offset] = deposit_loc { index };
      current_ptr[offset++]   = current;
    };

    store_current(
      i1,
      Vec3yee(
        Fx1 * (1.0 - Wy1) * (1.0 - Wz1),
        Fy1 * (1.0 - Wx1) * (1.0 - Wz1),
        Fz1 * (1.0 - Wx1) * (1.0 - Wy1)));

    store_current(
      i2,
      Vec3yee(
        Fx2 * (1.0 - Wy2) * (1.0 - Wz2),
        Fy2 * (1.0 - Wx2) * (1.0 - Wz2),
        Fz2 * (1.0 - Wx2) * (1.0 - Wy2)));

    store_current(
      i1 + Vec3uz(1, 0, 0),
      Vec3yee(0, Fy1 * Wx1 * (1.0 - Wz1), Fz1 * Wx1 * (1.0 - Wy1)));

    store_current(
      i2 + Vec3uz(1, 0, 0),
      Vec3yee(0, Fy2 * Wx2 * (1.0 - Wz2), Fz2 * Wx2 * (1.0 - Wy2)));

    store_current(
      i1 + Vec3uz(0, 1, 0),
      Vec3yee(Fx1 * Wy1 * (1.0 - Wz1), 0, Fz1 * (1.0 - Wx1) * Wy1));

    store_current(
      i2 + Vec3uz(0, 1, 0),
      Vec3yee(Fx2 * Wy2 * (1.0 - Wz2), 0, Fz2 * (1.0 - Wx2) * Wy2));

    store_current(
      i1 + Vec3uz(0, 0, 1),
      Vec3yee(Fx1 * (1.0 - Wy1) * Wz1, Fy1 * (1.0 - Wx1) * Wz1, 0));

    store_current(
      i2 + Vec3uz(0, 0, 1),
      Vec3yee(Fx2 * (1.0 - Wy2) * Wz2, Fy2 * (1.0 - Wx2) * Wz2, 0));

    store_current(i1 + Vec3uz(0, 1, 1), Vec3yee(Fx1 * Wy1 * Wz1, 0, 0));

    store_current(i2 + Vec3uz(0, 1, 1), Vec3yee(Fx2 * Wy2 * Wz2, 0, 0));

    store_current(i1 + Vec3uz(1, 0, 1), Vec3yee(0, Fy1 * Wx1 * Wz1, 0));

    store_current(i2 + Vec3uz(1, 0, 1), Vec3yee(0, Fy2 * Wx2 * Wz2, 0));

    store_current(i1 + Vec3uz(1, 1, 0), Vec3yee(0, 0, Fz1 * Wx1 * Wy1));

    store_current(i2 + Vec3uz(1, 1, 0), Vec3yee(0, 0, Fz2 * Wx2 * Wy2));
  });

  auto unique_deposit_locations = thrust::device_vector<std::array<std::size_t, 3>>(N);
  auto reduced_currents = thrust::device_vector<std::array<yee_value_type, 3>>(N);

  const auto uniq_dep_locs_ptr = thrust::make_transform_output_iterator(
    unique_deposit_locations.begin(),
    [](const deposit_loc& loc) { return loc.idx.data; });

  const auto reduced_currents_ptr = thrust::make_transform_output_iterator(
    reduced_currents.begin(),
    [](const Vec3yee& J) { return J.data; });


  w1.wait();

  thrust::sort_by_key(
    thrust::device,
    deposit_locations.begin(),
    deposit_locations.end(),
    currents.begin());

  const auto [key_end, _] = thrust::reduce_by_key(
    thrust::device,
    deposit_locations.begin(),
    deposit_locations.end(),
    currents.begin(),
    uniq_dep_locs_ptr,
    reduced_currents_ptr,
    thrust::equal_to<deposit_loc> {},
    thrust::plus<toolbox::Vec3<float>> {});

  const auto num_of_uniq_locs = key_end - uniq_dep_locs_ptr;

  unique_deposit_locations.resize(num_of_uniq_locs);
  reduced_currents.resize(num_of_uniq_locs);

  return emf2::YeeLattice::CurrentContributions { .locations =
                                                    std::move(unique_deposit_locations),
                                                  .currents =
                                                    std::move(reduced_currents) };
}

void
  ParticleContainer::current_zigzag_1st(
    runko::VecGrid<emf2::YeeLattice::value_type>& Jout,
    const std::array<value_type, 3> lattice_origo_coordinates,
    const value_type cfl) const
{
  using Vec3uz         = toolbox::Vec3<std::size_t>;
  using yee_value_type = emf2::YeeLattice::value_type;
  using Vec3yee        = toolbox::Vec3<yee_value_type>;
  using Vec3val        = toolbox::Vec3<value_type>;
  using Vec3d          = toolbox::Vec3<double>;

  const auto Jmds    = Jout.mds();
  const auto pos_mds = pos_.mds();
  const auto vel_mds = vel_.mds();
  const auto charge  = charge_;

  tyvi::mdgrid_work {}
    .for_each_index(
      pos_mds,
      [=](const auto idx) {
        // NOTE: performing velocity calculations via doubles to retain accuracy

        const auto u      = Vec3val(vel_mds[idx]).template as<double>();
        const auto invgam = 1.0 / std::sqrt(1.0 + toolbox::dot(u, u));

        const auto x2 = Vec3val(pos_mds[idx]).template as<double>() -
                        Vec3val(lattice_origo_coordinates).template as<double>();

        const auto x1 = x2 - cfl * invgam * u;

        // These are equivalent to floor assuming that x1 and x2 do not have negative
        // values.
        const auto i1 = x1.template as<std::size_t>();
        const auto i2 = x2.template as<std::size_t>();

        const auto relay = [&](const std::size_t j) -> double {
          const auto a  = static_cast<double>(min(i1(j), i2(j)) + 1);
          const auto b1 = static_cast<double>(max(i1(j), i2(j)));
          const auto b2 = 0.5 * (x1(j) + x2(j));
          const auto b  = max(b1, b2);
          return min(a, b);
        };

        const auto x_relay = Vec3d(relay(0), relay(1), relay(2));

        const auto F1 = charge * (x_relay - x1);
        const auto F2 = charge * (x2 - x_relay);

        const auto W1 = 0.5 * (x1 + x_relay) - i1.template as<double>();
        const auto W2 = 0.5 * (x2 + x_relay) - i2.template as<double>();

        const auto [Fx1, Fy1, Fz1] = F1.data;
        const auto [Fx2, Fy2, Fz2] = F2.data;
        const auto [Wx1, Wy1, Wz1] = W1.data;
        const auto [Wx2, Wy2, Wz2] = W2.data;

        const auto store_current = [&](const Vec3uz& index, const Vec3yee& current) {
          auto* const Jx = &thrust::raw_reference_cast(Jmds[index.data][0]);
          auto* const Jy = &thrust::raw_reference_cast(Jmds[index.data][1]);
          auto* const Jz = &thrust::raw_reference_cast(Jmds[index.data][2]);

          // See hip docs for these.
          std::ignore = ::unsafeAtomicAdd(Jx, current[0]);
          std::ignore = ::unsafeAtomicAdd(Jy, current[1]);
          std::ignore = ::unsafeAtomicAdd(Jz, current[2]);
        };

        store_current(
          i1,
          Vec3yee(
            Fx1 * (1.0 - Wy1) * (1.0 - Wz1),
            Fy1 * (1.0 - Wx1) * (1.0 - Wz1),
            Fz1 * (1.0 - Wx1) * (1.0 - Wy1)));

        store_current(
          i2,
          Vec3yee(
            Fx2 * (1.0 - Wy2) * (1.0 - Wz2),
            Fy2 * (1.0 - Wx2) * (1.0 - Wz2),
            Fz2 * (1.0 - Wx2) * (1.0 - Wy2)));

        store_current(
          i1 + Vec3uz(1, 0, 0),
          Vec3yee(0, Fy1 * Wx1 * (1.0 - Wz1), Fz1 * Wx1 * (1.0 - Wy1)));

        store_current(
          i2 + Vec3uz(1, 0, 0),
          Vec3yee(0, Fy2 * Wx2 * (1.0 - Wz2), Fz2 * Wx2 * (1.0 - Wy2)));

        store_current(
          i1 + Vec3uz(0, 1, 0),
          Vec3yee(Fx1 * Wy1 * (1.0 - Wz1), 0, Fz1 * (1.0 - Wx1) * Wy1));

        store_current(
          i2 + Vec3uz(0, 1, 0),
          Vec3yee(Fx2 * Wy2 * (1.0 - Wz2), 0, Fz2 * (1.0 - Wx2) * Wy2));

        store_current(
          i1 + Vec3uz(0, 0, 1),
          Vec3yee(Fx1 * (1.0 - Wy1) * Wz1, Fy1 * (1.0 - Wx1) * Wz1, 0));

        store_current(
          i2 + Vec3uz(0, 0, 1),
          Vec3yee(Fx2 * (1.0 - Wy2) * Wz2, Fy2 * (1.0 - Wx2) * Wz2, 0));

        store_current(i1 + Vec3uz(0, 1, 1), Vec3yee(Fx1 * Wy1 * Wz1, 0, 0));

        store_current(i2 + Vec3uz(0, 1, 1), Vec3yee(Fx2 * Wy2 * Wz2, 0, 0));

        store_current(i1 + Vec3uz(1, 0, 1), Vec3yee(0, Fy1 * Wx1 * Wz1, 0));

        store_current(i2 + Vec3uz(1, 0, 1), Vec3yee(0, Fy2 * Wx2 * Wz2, 0));

        store_current(i1 + Vec3uz(1, 1, 0), Vec3yee(0, 0, Fz1 * Wx1 * Wy1));

        store_current(i2 + Vec3uz(1, 1, 0), Vec3yee(0, 0, Fz2 * Wx2 * Wy2));
      })
    .wait();
}

}  // namespace pic2
