// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/pic/particle.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#include <array>
#include <boost/ut.hpp>
#include <cstddef>
#include <vector>

namespace {

using namespace boost::ut;
using pic::ParticleContainer;
using runko::dead_prtc_id;

constexpr auto X = ParticleContainer::no_source;

/// Per slot (id, x, vy, w) copied to the host; x, vy, w are read only for alive slots.
struct slots {
  thrust::host_vector<runko::prtc_id_type> id;
  thrust::host_vector<float> x, vy, w;
};

slots
  read(ParticleContainer& pc)
{
  const auto n   = pc.size();
  const auto ids = pc.ids_mds();
  const auto pos = pc.pos_mds();
  const auto vel = pc.vel_mds();
  const auto wgt = pc.wgt_mds();
  const auto weighted = pc.has_particle_weights();
  const auto i0  = thrust::counting_iterator<std::size_t>(0uz);
  const auto i1  = thrust::counting_iterator<std::size_t>(n);
  auto id        = thrust::device_vector<runko::prtc_id_type>(n);
  auto x         = thrust::device_vector<float>(n);
  auto vy        = thrust::device_vector<float>(n);
  auto w         = thrust::device_vector<float>(n);
  thrust::transform(i0, i1, id.begin(), [=](std::size_t i) { return ids[i][]; });
  thrust::transform(i0, i1, x.begin(), [=](std::size_t i) { return pos[i][0]; });
  thrust::transform(i0, i1, vy.begin(), [=](std::size_t i) { return vel[i][1]; });
  thrust::transform(i0, i1, w.begin(), [=](std::size_t i) {
    return weighted ? wgt[i][] : 1.0f;
  });
  return { id, x, vy, w };
}

/// Slot j must hold the particle k = want[j] (id 10 + k, x = k, vy = 2k, w = 3k + 1),
/// or be dead when want[j] < 0.
void
  expect_slots(ParticleContainer& pc, const std::vector<int>& want, const bool weighted)
{
  const auto s = read(pc);
  expect(that % pc.size() == want.size());
  for(auto j = 0uz; j < want.size(); ++j) {
    const auto k = want[j];
    if(k < 0) {
      expect(that % s.id[j] == dead_prtc_id);
      continue;
    }
    expect(that % s.id[j] == runko::prtc_id_type(10 + k));
    expect(that % s.x[j] == float(k));
    expect(that % s.vy[j] == float(2 * k));
    if(weighted) { expect(that % s.w[j] == float(3 * k + 1)); }
  }
}

[[maybe_unused]]
const suite<"particle container"> _ = [] {
  "reorder permutes, grows with dead gaps and shrinks; append and set"_test = [] {
    for(const auto mass: { 1.0, 0.0 }) {  // a massless species carries weights
      auto pc = ParticleContainer({ .N = 0, .charge = 0, .mass = mass });
      auto p  = std::vector<runko::ParticleState<double>> {};
      for(auto k = 0; k < 4; ++k) {
        p.push_back({ .pos { double(k), 0, 0 },
                      .vel { 0, double(2 * k), 0 },
                      .id  = runko::prtc_id_type(10 + k),
                      .wgt = double(3 * k + 1) });
      }
      pc.add_particles(p);
      const auto weighted = pc.has_particle_weights();
      const auto n0       = pc.reorder_count();

      pc.reorder(thrust::device_vector<runko::index_t>(std::vector { 3u, 1u, 0u, 2u }));
      expect_slots(pc, { 3, 1, 0, 2 }, weighted);

      pc.reorder(thrust::device_vector<runko::index_t>(std::vector { 0u, X, 1u, X, 2u, 3u }));
      expect_slots(pc, { 3, -1, 1, -1, 0, 2 }, weighted);
      expect(that % pc.number_alive() == 4uz);

      // a dead slot moves as dead; particle 2 is dropped
      pc.reorder(thrust::device_vector<runko::index_t>(std::vector { 4u, 2u, 0u, 1u }));
      expect_slots(pc, { 0, 1, 3, -1 }, weighted);
      expect(that % pc.number_alive() == 3uz);
      expect(that % pc.reorder_count() == n0 + 3);

      // append fills the dead slots behind the last alive one, and the container ends there
      pc.add_particles(std::vector { p[2], p[2] });
      expect_slots(pc, { 0, 1, 3, 2, 2 }, weighted);
      expect(that % pc.reorder_count() == n0 + 3);

      // set replaces every slot: a reorder
      const auto other = ParticleContainer(pc);
      pc.set(std::array { ParticleContainer::specific_span { .span      = { 1, 3 },
                                                             .container = &other } });
      expect_slots(pc, { 1, 3 }, weighted);
      expect(that % pc.reorder_count() == n0 + 4);
    }
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{ return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv })); }
