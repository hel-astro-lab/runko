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

/// Per slot (id, x, vy, payload a, b) copied to the host; all but id are read only for
/// alive slots.
struct slots {
  thrust::host_vector<runko::prtc_id_type> id;
  thrust::host_vector<float> x, vy, a, b;
};

slots
  read(ParticleContainer& pc)
{
  const auto n   = pc.size();
  const auto ids = pc.ids_mds();
  const auto pos = pc.pos_mds();
  const auto vel = pc.vel_mds();
  const auto i0  = thrust::counting_iterator<std::size_t>(0uz);
  const auto i1  = thrust::counting_iterator<std::size_t>(n);
  auto id        = thrust::device_vector<runko::prtc_id_type>(n);
  auto x         = thrust::device_vector<float>(n);
  auto vy        = thrust::device_vector<float>(n);
  auto a         = thrust::device_vector<float>(n);
  auto b         = thrust::device_vector<float>(n);
  const auto pa  = pc.payload_mds(0);
  const auto pb  = pc.payload_mds(1);
  thrust::transform(i0, i1, a.begin(), [=](std::size_t i) { return pa[i][]; });
  thrust::transform(i0, i1, b.begin(), [=](std::size_t i) { return pb[i][]; });
  thrust::transform(i0, i1, id.begin(), [=](std::size_t i) { return ids[i][]; });
  thrust::transform(i0, i1, x.begin(), [=](std::size_t i) { return pos[i][0]; });
  thrust::transform(i0, i1, vy.begin(), [=](std::size_t i) { return vel[i][1]; });
  return { id, x, vy, a, b };
}

/// Slot j must hold the particle k = want[j] (id 10 + k, x = k, vy = 2k, payload
/// (4k + 1, -k)), or be dead when want[j] < 0.
void
  expect_slots(ParticleContainer& pc, const std::vector<int>& want)
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
    expect(that % s.a[j] == float(4 * k + 1));
    expect(that % s.b[j] == float(-k));
  }
}

[[maybe_unused]]
const suite<"particle container"> _ = [] {
  "reorder permutes, grows with gaps and shrinks; append, set; payload"_test = [] {
    auto pc = ParticleContainer({ .N = 0, .charge = 0, .mass = 1, .n_payload = 2 });
    auto p  = std::vector<runko::ParticleState<double>> {};
    auto pl = std::vector<float> {};  // two components per particle
    for(auto k = 0; k < 4; ++k) {
      p.push_back({ .pos { double(k), 0, 0 },
                    .vel { 0, double(2 * k), 0 },
                    .id = runko::prtc_id_type(10 + k) });
      pl.insert(pl.end(), { float(4 * k + 1), float(-k) });
    }
    pc.add_particles(p, pl);
    const auto n0      = pc.reorder_count();
    const auto reorder = [&](const std::vector<runko::index_t>& src) {
      pc.reorder(thrust::device_vector<runko::index_t>(src));
    };

    reorder({ 3, 1, 0, 2 });
    expect_slots(pc, { 3, 1, 0, 2 });

    reorder({ 0, X, 1, X, 2, 3 });
    expect_slots(pc, { 3, -1, 1, -1, 0, 2 });
    expect(that % pc.number_alive() == 4uz);

    // a dead slot moves as dead; particle 2 is dropped
    reorder({ 4, 2, 0, 1 });
    expect_slots(pc, { 0, 1, 3, -1 });
    expect(that % pc.number_alive() == 3uz);
    expect(that % pc.reorder_count() == n0 + 3);

    // append fills the dead slots behind the last alive one, and the container ends
    // there
    pc.add_particles(std::vector { p[2], p[2] }, std::vector { 9.f, -2.f, 9.f, -2.f });
    expect_slots(pc, { 0, 1, 3, 2, 2 });
    expect(that % pc.reorder_count() == n0 + 3);

    // append_dead keeps the particles and their payload
    pc.append_dead(3);
    expect_slots(pc, { 0, 1, 3, 2, 2, -1, -1, -1 });

    // a payload of the wrong size, or none, is refused
    const auto one = std::vector { p[0] };
    expect(throws([&] { pc.add_particles(one); }));
    expect(throws([&] { pc.add_particles(one, std::vector { 1.f }); }));

    // set replaces every slot: a reorder
    auto other = ParticleContainer(pc);
    pc.set(std::array { ParticleContainer::specific_span { .span      = { 1, 3 },
                                                           .container = &other } });
    expect_slots(pc, { 1, 3 });
    expect(that % pc.reorder_count() == n0 + 4);
    expect_slots(other, { 0, 1, 3, 2, 2, -1, -1, -1 });  // a copy

    auto assigned = ParticleContainer({ .N = 0, .charge = 0, .mass = 1 });
    assigned      = other;
    expect_slots(assigned, { 0, 1, 3, 2, 2, -1, -1, -1 });
  };
};

}  // namespace

int
  main(int argc, const char** argv)
{ return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv })); }
