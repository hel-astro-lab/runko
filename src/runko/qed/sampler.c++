// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

// The step of the sampler: layout cadence, channel order, systematic rounding, the
// placement of products, statistics.

#include "runko/qed/sampler.h"

#include "runko/qed/interactions/qed1_bw.h"
#include "runko/qed/interactions/qed1_syn.h"
#include "runko/qed/interactions/qed2_ann.h"
#include "runko/qed/interactions/qed2_comp.h"
#include "runko/qed/interactions/qed2_gg.h"
#include "runko/qed/interactions/qed2_ssa.h"
#include "runko/qed/weights.h"
#include "runko/tools/math.h"
#include "runko/tools/random.h"
#include "thrust/binary_search.h"
#include "thrust/copy.h"
#include "thrust/fill.h"
#include "thrust/for_each.h"
#include "thrust/host_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/reduce.h"
#include "thrust/scan.h"
#include "thrust/sort.h"
#include "thrust/transform.h"
#include "thrust/transform_reduce.h"
#include "tyvi/mdgrid.h"

#include <algorithm>
#include <format>
#include <sstream>
#include <stdexcept>

namespace qed {

namespace {

// in Channel order; the keys of counts() and moved()
constexpr const char* NAMES[] = { "compton-", "compton+", "pairann", "photann", "ssa-",
                                  "ssa+",     "syn-",     "syn+",    "bw" };

}  // namespace

const char*
  channel_name(const Channel c)
{ return NAMES[static_cast<int>(c)]; }

bool
  is_twobody(const Channel c)
{ return static_cast<int>(c) <= static_cast<int>(Channel::ssa_p); }

Sampler::Sampler(
  const Box& b,
  const Params& p,
  const std::vector<std::string>& channels) :
  box { b },
  par { p }
{
  // channels from ssa_m on need a field, b_perp > 0; a given list restricts further
  const auto field = b.b_perp() > 0.0;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    const auto c    = static_cast<Channel>(i);
    const auto need = static_cast<int>(c) >= static_cast<int>(Channel::ssa_m);
    auto on         = not need or field;
    if(not channels.empty()) {
      on = on and std::ranges::find(channels, std::string(NAMES[i])) != channels.end();
    }
    enabled[i] = on;
  }

  nstaged_.assign(1, 0);
  placed_.assign(1, Placed {});
  dstat_.assign(N_CHANNELS, DeviceStats {});
  tally_.assign(1, Tally {});
}

thrust::device_vector<double>&
  Sampler::carry(const int l, const int role, const int ncell)
{
  auto& car = round_.carry[3 * static_cast<std::size_t>(l) + static_cast<std::size_t>(role)];

  if(car.size() != static_cast<std::size_t>(ncell)) {
    auto h   = thrust::host_vector<double>(static_cast<std::size_t>(ncell));
    auto rng = runko::Rng { runko::rng_key(par.seed, 0, 0, stream::carry(l, role)), 0 };
    for(auto& x: h) { x = rng.uniform(); }
    car = h;
  }

  return car;
}

// Systematic rounding: each event is decided against a running remainder of its cell,
// not an independent uniform; the marginals are unchanged and a cell's summed error
// stays below one. In cell order a cell's counts telescope to floor(base + sum) -
// floor(base), and its new remainder is frac(base + sum).
void
  Sampler::select(
    const tyvi::mdgrid_work& w,
    bool sorted,
    const int l,
    const int role,
    const int ncell,
    const std::size_t n,
    const std::int32_t* cell,
    const double* q,
    std::int32_t* out)
{
  if(n == 0) { return; }

  auto& car        = carry(l, role, ncell);
  const auto car_p = thrust::raw_pointer_cast(car.data());  // remainder per cell

  // one buffer, of the pass's event count, serves every select of the pass: no resize
  // may reallocate under a select still in flight on the queue
  if(round_.cs.size() < n) { round_.cs.resize(n); }
  const auto cs   = ptr(round_.cs);
  const auto ord0 = thrust::counting_iterator<I>(0);
  const auto nn   = static_cast<std::ptrdiff_t>(n);
  sorted          = sorted or ncell == 1;  // one cell: every order is cell-sorted

  // the events in cell order: perm (null for the identity) and their cells along it
  const I* perm            = nullptr;
  const std::int32_t* keys = cell;
  if(not sorted) {
    round_.perm.resize(n);
    round_.keys.resize(n);
    thrust::copy(w.on_this(), cell, cell + nn, round_.keys.begin());
    thrust::copy(w.on_this(), ord0, ord0 + nn, round_.perm.begin());
    thrust::stable_sort_by_key(
      w.on_this(),
      round_.keys.begin(),
      round_.keys.end(),
      round_.perm.begin());
    perm = ptr(round_.perm);
    keys = ptr(round_.keys);
  }

  // event at sorted position i
  const auto at = [=](const I i) { return perm ? perm[i] : i; };
  const auto qs =
    thrust::make_transform_iterator(ord0, [=](const I i) { return q[at(i)]; });
  thrust::inclusive_scan_by_key(w.on_this(), keys, keys + nn, qs, cs);

  thrust::for_each(w.on_this(), ord0, ord0 + nn, [=](const I i) {
    const auto j    = at(i);
    const auto base = car_p[keys[i]];
    out[j]          = static_cast<std::int32_t>(
      sstd::floor(base + cs[i]) - sstd::floor(base + cs[i] - q[j]));
  });

  // a second pass, since every event of the cell reads the old remainder above: the
  // cell's last event advances it
  thrust::for_each(w.on_this(), ord0, ord0 + nn, [=](const I i) {
    if(static_cast<std::size_t>(i) + 1 < n and keys[i + 1] == keys[i]) { return; }
    const auto x   = car_p[keys[i]] + cs[i];
    car_p[keys[i]] = x - sstd::floor(x);
  });
}

// ---------------------------------------------------------------- the step
void
  Sampler::layout(TileView& t)
{
  const auto lp = LayoutParams { .nbin       = par.nbin,
                                 .ncell_edge = par.ncell_edge,
                                 .slack      = par.slack,
                                 .smin       = par.smin };

  for(auto s = 0uz; s < 3uz; ++s) {
    t.layout[s]
      ->build(*t.con[s], lp, t.mins, t.maxs, runko::rng_key(par.seed, istep_, s, stream::layout));
  }
}

// a species loses its tables when a source fills a fresh or reset tile
void
  Sampler::lay_if_needed(TileView& t, const bool epoch)
{
  auto need = epoch;
  for(auto s = 0uz; s < 3uz; ++s) { need = need or not t.layout[s]->laid(*t.con[s]); }
  if(need) { this->layout(t); }
}

void
  Sampler::step(TileView& t)
{
  // ---- layout every nsort steps, or as soon as one went stale
  const auto epoch = istep_ % static_cast<std::uint64_t>(par.nsort) == 0;
  this->lay_if_needed(t, epoch);

  if(not epoch) {
    // the PIC push moved particles without widening the key bounds; this sweep keeps
    // every envelope a majorant
    const auto w = tyvi::mdgrid_work {};
    for(auto s = 0uz; s < 3uz; ++s) { t.layout[s]->touch_all(*t.con[s], w); }
    w.wait();
  }

  // ---- context and slot statistics
  Ctx c { t, par.nbin };
  for(auto s = 0uz; s < 3uz; ++s) {
    stats.slots[s] += static_cast<double>(t.con[s]->size());
    stats.live[s] += static_cast<double>(t.con[s]->number_alive());
    t.ids_used[s] = 0;
  }

  // ---- the channels in random order, each reading the states the previous one left;
  // the splitting error is second order in the rates and cancels in expectation
  auto order = std::vector<int> {};
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { order.push_back(static_cast<int>(i)); }
  }

  // Fisher-Yates shuffle
  auto rng = runko::Rng { runko::rng_key(par.seed, istep_, 0, stream::channel_order), 0 };
  for(auto i = order.size(); i > 1; --i) {
    const auto j = static_cast<std::size_t>(rng.integer(i));
    std::swap(order[i - 1], order[j]);
  }

  for(const auto l: order) {
    const auto which = static_cast<Channel>(l);
    if(is_twobody(which)) {
      twobody(c, which);
    } else if(which == Channel::bw) {
      breit_wheeler(c);
    } else {
      synchrotron(c, which);
    }
  }

  ++stats.steps;
  ++istep_;
}

// ---------------------------------------------------------------- products
Stager
  Sampler::stager(const std::size_t n)
{
  staged_.resize(n);
  nstaged_[0] = 0;
  return Stager { ptr(staged_), ptr(nstaged_) };
}

void
  Sampler::place(Ctx& c)
{
  const auto nst = static_cast<std::size_t>(nstaged_[0]);  // staged records
  if(nst == 0) { return; }

  const auto w   = tyvi::mdgrid_work {};
  const auto stg = thrust::raw_pointer_cast(staged_.data());
  const auto r0  = thrust::counting_iterator<I>(0);

  // copies per species, then room for them at the end of every container
  std::array<std::size_t, 3> tot {};
  for(auto s = 0; s < 3; ++s) {
    tot[static_cast<std::size_t>(s)] =
      static_cast<std::size_t>(thrust::transform_reduce(
        w.on_this(),
        r0,
        r0 + static_cast<std::ptrdiff_t>(nst),
        [=](const I k) {
          return stg[k].species == s ? static_cast<long>(stg[k].n_r) : 0L;
        },
        0L,
        thrust::plus<long> {}));
  }

  // the tail of every species starts at its old container end
  std::array<I, 3> tail_base {};
  for(auto s = 0uz; s < 3uz; ++s) {
    tail_base[s] = static_cast<I>(c.t->con[s]->size());
    if(tot[s] > 0) {
      c.t->con[s]->append_dead(tot[s]);
      c.t->layout[s]->resized(*c.t->con[s]);
    }
  }
  c.refresh_all();  // the containers were resized

  placed_[0]     = Placed {};
  const auto cnt = ptr(placed_);
  const auto v   = c.v;
  const auto idb = c.t->id_base;  // first free id per species

  // one lane per staged record; each copy takes a slack slot of its key or a tail slot
  thrust::for_each(
    w.on_this(),
    r0,
    r0 + static_cast<std::ptrdiff_t>(nst),
    [=](const I k) {
      const auto rec = stg[k];
      const auto s   = static_cast<std::size_t>(rec.species);
      const auto pv  = v[s];
      const auto lay = pv.lay;
      const auto e   = energy_of(pv.mass, rec.u[0], rec.u[1], rec.u[2]);
      const auto kk  = lay.key_of(rec.cell, lay.bin_of(e), rec.bit);

      for(auto n = 0; n < rec.n_r; ++n) {
        const auto id =
          idb[s] +
          static_cast<runko::prtc_id_type>(
            sstd::atomic_fetch_add(&cnt->id[rec.species], std::int32_t { 1 }));

        // the copy's rank behind what the key already holds
        const auto f    = sstd::atomic_fetch_add(&lay.fill[kk], I { 1 });
        const auto at   = lay.b[kk] + f;
        const auto fits = f < lay.cap[kk];
        I slot;

        // claim the slot by CAS on its dead id: an arrival of the particle exchange
        // may have taken it since the layout
        if(
          fits and
          sstd::atomic_compare_exchange(&pv.ids[at][], runko::dead_prtc_id, id)) {
          slot = at;  // reactive at once

          // widen the key's bounds, so the envelopes of the channels still to come this
          // step stay upper bounds
          lay.widen(kk, e, rec.w);
          sstd::atomic_add(&cnt->slack, std::int32_t { 1 });
        } else {
          if(fits) {
            // the arrival stays and the fill now counts it in the key: widen by it too
            lay.widen(kk, static_cast<float>(pv.e(at)), static_cast<float>(pv.w(at)));
            sstd::atomic_add(&cnt->held, std::int32_t { 1 });
          }
          const auto t =
            sstd::atomic_fetch_add(&cnt->tail[rec.species], std::int32_t { 1 });
          slot = tail_base[s] + static_cast<I>(t);  // inert until the next sort
          sstd::atomic_add(&cnt->tailn, std::int32_t { 1 });
        }

        sstd::atomic_add(&lay.prod[kk], I { 1 });  // demand, whether or not it fitted

        for(auto q = 0uz; q < 3uz; ++q) {
          pv.pos[slot][q] = rec.pos[q];
          pv.vel[slot][q] = rec.u[q];
        }
        pv.ids[slot][] = id;
        if(pv.weighted) { pv.wgt[slot][] = rec.w; }
      }
    });
  w.wait();

  {
    // ids consumed and placement statistics
    const Placed h = placed_[0];
    for(auto s = 0uz; s < 3uz; ++s) {
      c.t->ids_used[s] += static_cast<std::size_t>(h.id[s]);
      c.t->id_base[s] += static_cast<runko::prtc_id_type>(h.id[s]);
    }
    stats.slackprod += h.slack;
    stats.tailprod += h.tailn;
    stats.heldprod += h.held;
  }
}

// Each field is written by one side only, so the sum is exact; fired is the one-body
// accepted count.
ChannelStats
  Sampler::channel_snap(
    const std::size_t l,
    const thrust::host_vector<DeviceStats>& h) const
{
  auto s        = stats.ch[l];
  const auto& d = h[l];
  s.acc += d.fired;
  s.facc += d.facc;
  s.wem += d.wem;
  s.clamped += d.clamped;
  s.shbar_num += d.shnum;
  s.shbar_den += d.shden;
  return s;
}

// accepted events per enabled channel, cumulative
std::map<std::string, double>
  Sampler::counts() const
{
  const thrust::host_vector<DeviceStats> h = dstat_;
  std::map<std::string, double> m;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { m[NAMES[i]] = channel_snap(i, h).acc; }
  }
  return m;
}

// weight moved per enabled channel, cumulative: sum f min(w_i, w_j) over the events
std::map<std::string, double>
  Sampler::moved() const
{
  const thrust::host_vector<DeviceStats> h = dstat_;
  std::map<std::string, double> m;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(enabled[i]) { m[NAMES[i]] = channel_snap(i, h).wmoved; }
  }
  return m;
}

std::string
  Sampler::report() const
{
  std::ostringstream o;
  const auto ns =
    static_cast<double>(std::max<std::size_t>(1, stats.steps));  // everything per step

  o << std::format(
    "kernels per step over {} steps ({} bin(s), sort every {}, slack {}):\n",
    stats.steps,
    par.nbin,
    par.nsort,
    par.slack);
  o << std::format(
    "   slots / live LPs       : ph={:.0f}/{:.0f}, e-={:.0f}/{:.0f}, "
    "e+={:.0f}/{:.0f}\n",
    stats.slots[0] / ns,
    stats.live[0] / ns,
    stats.slots[1] / ns,
    stats.live[1] / ns,
    stats.slots[2] / ns,
    stats.live[2] / ns);
  o << std::format("   cell-kernel lanes      : {:.0f}\n", stats.cell_lanes / ns);

  const thrust::host_vector<DeviceStats> h = dstat_;
  for(auto i = 0uz; i < N_CHANNELS; ++i) {
    if(not enabled[i]) { continue; }
    const auto s = channel_snap(i, h);

    if(is_twobody(static_cast<Channel>(i))) {
      o << std::format(
        "   {:9s}: lanes {:.1f}, virtual {:.1f}%, weight-stage {:.3f}, accepted "
        "{:.2f}, <sigmahat> {:.4f}, <f> {:.2f}, nsub {} ({} steps), clamped {:.0f}\n",
        NAMES[i],
        s.M / ns,
        100.0 * s.virt / std::max(1.0, s.M),
        s.wpass / std::max(1.0, s.M - s.virt),
        s.acc / ns,
        s.shbar_num / std::max(1e-30, s.shbar_den),
        s.facc / std::max(1.0, s.acc),
        s.nsub_max_seen,
        s.nsub_steps,
        s.clamped);
    } else {
      // one lane per slot: the dead share is 1 - live / slots of the channel's species
      const auto is   = static_cast<Channel>(i) == Channel::syn_m   ? 1uz
                        : static_cast<Channel>(i) == Channel::syn_p ? 2uz
                                                                    : 0uz;
      const auto dead = 100.0 * (1.0 - stats.live[is] / std::max(stats.slots[is], 1.0));
      o << std::format(
        "   {:9s}: lanes {:.0f}, events {:.2f} (dead {:.1f}%), emitted weight/step "
        "{:.3g}, <f> {:.2f}\n",
        NAMES[i],
        s.M / ns,
        s.acc / ns,
        dead,
        s.wem / ns,
        s.facc / std::max(1.0, s.acc));
    }
  }

  o << std::format(
    "   products slack / tail  : {:.2f}/{:.2f} ({:.2f} found an arrival in slack)\n",
    stats.slackprod / ns,
    stats.tailprod / ns,
    stats.heldprod / ns);
  return o.str();
}

}  // namespace qed
