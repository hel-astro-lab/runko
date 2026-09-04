// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/qed/layout.h"

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
#include "tyvi/mdgrid.h"

#include <limits>
#include <stdexcept>

namespace qed {

using I = runko::index_t;

Pop::Pop(const bool keyed, const double mass) : keyed_ { keyed }, mass_ { static_cast<float>(mass) } {}

void
  Pop::layout(
    pic::ParticleContainer& con,
    const LayoutParams& par,
    const std::array<float, 3> mins,
    const std::array<float, 3> maxs,
    const std::uint64_t rng_key)
{
  // The only kernel of the scheme that runs over the particles; it is amortized over the sort
  // epoch and shared with the sort the particle pusher needs anyway. The bins are log-spaced
  // in the kinetic energy between the extremes of the live population, so every bound is
  // equally tight in every bin, and a sparse tail costs nothing: an empty bin launches no lane.
  par_  = par;  // the bin count, cell count and slack rule of this epoch
  mins_ = mins;
  maxs_ = maxs;
  // reaction sub-cells of the tile
  const auto ncell = par.ncell_edge * par.ncell_edge * par.ncell_edge;
  // keys: cells x bins, x 2 photon bits
  nkeys_           = (keyed_ ? 2 : 1) * ncell * par.nbin;

  const auto n = con.size();
  const auto w = tyvi::mdgrid_work {};

  // ---- device views of the container and the two per-slot predicates the kernels share
  const auto ids_mds = con.ids_mds();
  const auto vel_mds = con.vel_mds();
  const auto pos_mds = con.pos_mds();
  const auto m       = mass_;  // rest mass of the species: 0 for photons

  // energy of slot i: |k| for photons, gamma for leptons
  const auto energy = [=](const I i) {
    const auto vx = vel_mds[i][0], vy = vel_mds[i][1], vz = vel_mds[i][2];
    return sstd::sqrt(m * m + vx * vx + vy * vy + vz * vz);
  };
  // a slot holds a live LP unless its id is the dead marker
  const auto alive = [=](const I i) { return ids_mds[i][] != runko::dead_prtc_id; };

  const auto ord0 = thrust::counting_iterator<I>(0);  // slot indices 0 .. n-1
  const auto ordN = ord0 + static_cast<std::ptrdiff_t>(n);

  // ---- energy extremes of the live population fix the bins of the epoch
  // sentinel: dead slots take no part in the extremes
  const auto inf  = std::numeric_limits<float>::infinity();
  // kinetic energy of a live slot, +inf (min) or -inf (max) for a dead one
  const auto kin  = [=](const I i) {
    return alive(i) ? energy(i) - m : inf;
  };
  const auto kin_hi = [=](const I i) {
    return alive(i) ? energy(i) - m : -inf;
  };
  // smallest and largest kinetic energy of the live population (two reductions)
  const auto lo = thrust::reduce(
    w.on_this(),
    thrust::make_transform_iterator(ord0, kin),
    thrust::make_transform_iterator(ordN, kin),
    inf,
    thrust::minimum<float> {});
  const auto hi = thrust::reduce(
    w.on_this(),
    thrust::make_transform_iterator(ord0, kin_hi),
    thrust::make_transform_iterator(ordN, kin_hi),
    -inf,
    thrust::maximum<float> {});
  // log-spaced bin edges between the extremes, widened so the maximum falls inside the last bin
  if(par.nbin > 1 and lo < hi) {
      // floor against a lepton exactly at rest
    const auto vlo = sstd::max(static_cast<double>(lo), 1e-30);
    auto vhi       = static_cast<double>(hi) * (1.0 + 1e-6);
    vhi            = sstd::max(vhi, vlo * (1.0 + 1e-6));
    llo_           = sstd::log(vlo);  // log of the lower edge
    dl_            = (sstd::log(vhi) - llo_) / par.nbin;  // log width of one bin
  } else {
    llo_ = -1e300;  // everything in bin 0
    dl_  = 1e300;
  }

  // ---- key of every slot, dead ones to a sentinel past the last key
  // a device view of the tables; only its bin/cell/key rules are used here
  auto v = this->view();
  v.llo  = llo_;
  v.dl   = dl_;
  keys_.resize(n);  // key per slot
  ord_.resize(n);  // slot index per key, the permutation the sort produces
  const auto sentinel = static_cast<I>(nkeys_);  // dead slots sort past every key
  const auto keyed    = keyed_;
  // one lane per slot: its key from the reaction cell, the energy bin and a fresh photon bit
  thrust::transform(w.on_this(), ord0, ordN, keys_.begin(), [=](const I i) -> I {
    if(not alive(i)) { return sentinel; }
    const auto cell = v.cell_of(pos_mds[i][0], pos_mds[i][1], pos_mds[i][2]);
    const auto ebin = v.bin_of(static_cast<double>(energy(i)));
    const auto bit  = keyed ? static_cast<int>(runko::Rng { rng_key, i }.uniform() < 0.5) : 0;
    return static_cast<I>(v.key_of(cell, ebin, bit));
  });
  // sort the slot indices by key: ord_ then lists the live slots key by key
  thrust::copy(w.on_this(), ord0, ordN, ord_.begin());
  thrust::sort_by_key(w.on_this(), keys_.begin(), keys_.end(), ord_.begin());

  // ---- counts per key from the sorted keys, then the capacities and offsets
  const auto nk = static_cast<std::size_t>(nkeys_);  // number of keys
  cnt_.resize(nk + 1);  // starts of every key in the sorted list, plus the end
  const auto key0 = thrust::counting_iterator<I>(0);
  // first sorted position of every key (a vectorized binary search over the sorted keys)
  thrust::lower_bound(
    w.on_this(),
    keys_.begin(),
    keys_.end(),
    key0,
    key0 + static_cast<std::ptrdiff_t>(nk + 1),
    cnt_.begin());  // cnt_[k] = first sorted position of key k; cnt_[nk] = live total
  // starts[k+1] - starts[k] = live LPs of key k
  const auto starts = thrust::raw_pointer_cast(cnt_.data());
  // the product demand of the previous epoch, zero on the first layout
  if(prod_.size() != nk) {
    prod_.assign(nk, I { 0 });
  }
  cap_.resize(nk);  // slots of every key: live + slack
  b_.resize(nk);  // first slot of every key
  fill_.resize(nk);  // live LPs of every key, later live + appended
  const auto prod_p = thrust::raw_pointer_cast(prod_.data());
  const auto slack  = par.slack;  // slack as a fraction of the live LPs
  const auto smin   = static_cast<I>(par.smin);  // and its floor
  // one lane per key: its capacity
  thrust::transform(
    w.on_this(),
    key0,
    key0 + static_cast<std::ptrdiff_t>(nk),
    cap_.begin(),
    [=](const I k) -> I {
      // slack: a fraction of the live LPs, a floor, and what the key took in the previous
      // epoch with room to spare, so that a key that fills fast gets the slots it will need
      const auto c    = starts[k + 1] - starts[k];
      const auto need = static_cast<I>(sstd::ceil(1.5 * static_cast<double>(prod_p[k])));
      const auto sl   = static_cast<I>(sstd::ceil(slack * static_cast<double>(c)));
      return c + sstd::max(sstd::max(sl, smin), need);
    });
  // offsets from the capacities; the total is the padded size of the container
  thrust::exclusive_scan(w.on_this(), cap_.begin(), cap_.end(), b_.begin(), I { 0 });
  const auto total = thrust::reduce(w.on_this(), cap_.begin(), cap_.end(), I { 0 });
  // one lane per key: its live count starts the fill
  thrust::transform(
    w.on_this(),
    key0,
    key0 + static_cast<std::ptrdiff_t>(nk),
    fill_.begin(),
    [=](const I k) -> I { return starts[k + 1] - starts[k]; });
  // the demand of the new epoch starts from zero
  thrust::fill(w.on_this(), prod_.begin(), prod_.end(), I { 0 });

  // ---- slot of every live LP: its key's offset plus its rank inside the key
  slots_.resize(n);  // target slot of every current slot (dead ones unused)
  const auto keys_p  = thrust::raw_pointer_cast(keys_.data());
  const auto ord_p   = thrust::raw_pointer_cast(ord_.data());
  const auto b_p     = thrust::raw_pointer_cast(b_.data());
  const auto slots_p = thrust::raw_pointer_cast(slots_.data());
  // live LPs, the head of the sorted list
  const auto nlive   = static_cast<std::ptrdiff_t>(cnt_.back());
  // one lane per live LP in sorted order: rank j - starts[k] inside its key, after the key's offset
  thrust::for_each(w.on_this(), ord0, ord0 + nlive, [=](const I j) {
    const auto k       = keys_p[j];
    slots_p[ord_p[j]]  = b_p[k] + (j - starts[k]);
  });
  // scatter the container into the padded layout; every other slot becomes dead slack
  // (the wait: place() reads slots_ on its own streams, unordered against w)
  w.wait();
  con.place(slots_, static_cast<std::size_t>(total));

  // ---- key of every slot, and the bounds of every key from what it holds
  // key of every slot of the padded layout
  skey_.resize(static_cast<std::size_t>(total));
  const auto slot0 = thrust::counting_iterator<I>(0);
  // key of a slot = number of key offsets <= slot, minus one (vectorized binary search)
  thrust::upper_bound(
    w.on_this(),
    b_.begin(),
    b_.end(),
    slot0,
    slot0 + static_cast<std::ptrdiff_t>(total),
    skey_.begin());
  // upper_bound counts offsets <= slot; the key index is one less
  thrust::transform(
    w.on_this(),
    skey_.begin(),
    skey_.end(),
    skey_.begin(),
    [](const std::int32_t k) { return k - 1; });

  // bounds start empty and are widened by every live LP in touch_all
  wmax_.assign(nk, 0.0f);
  emin_.assign(nk, inf);
  emax_.assign(nk, 0.0f);
  laid_ = true;  // the tables are valid for the current container
  w.wait();      // touch_all reads skey_ on its own stream: the writers must be done
  this->touch_all(con);
}

void
  Pop::touch_all(pic::ParticleContainer& con)
{
  if(not laid_) { return; }  // no tables to widen; the next layout builds them
  // slots appended since the layout are tail (key -1) and widen nothing
  this->follow_container(con);
  const auto ids_mds = con.ids_mds();
  const auto vel_mds = con.vel_mds();
  const auto wgt_mds = con.wgt_mds();
  // ---- views and constants the bounds kernel needs
  const auto weighted = con.has_particle_weights();
  const auto w0       = static_cast<float>(con.weight());  // weight of a lepton species
  const auto m        = mass_;
  const auto v        = this->view();
  // one lane per slot: a live LP with a key widens the weight and energy bounds of that key
  tyvi::mdgrid_work {}
    .for_each_index(
      ids_mds,
      [=](const auto idx) {
        if(ids_mds[idx][] == runko::dead_prtc_id) { return; }
        const auto k = v.skey[idx[0]];
        if(k < 0) { return; }  // tail slot: no key yet
        const auto vx = vel_mds[idx][0], vy = vel_mds[idx][1], vz = vel_mds[idx][2];
        const auto e  = sstd::sqrt(m * m + vx * vx + vy * vy + vz * vz);
        const auto wg = weighted ? wgt_mds[idx][] : w0;  // its weight
        sstd::atomic_max(&v.wmax[k], wg);
        sstd::atomic_max(&v.emax[k], e);
        sstd::atomic_min(&v.emin[k], e);
      })
    .wait();
}

LayoutView
  Pop::view()
{
  // raw device pointers of the tables; valid until the next resize
  const auto ncell = par_.ncell_edge * par_.ncell_edge * par_.ncell_edge;
  return LayoutView { .skey       = thrust::raw_pointer_cast(skey_.data()),
                      .b          = thrust::raw_pointer_cast(b_.data()),
                      .cap        = thrust::raw_pointer_cast(cap_.data()),
                      .fill       = thrust::raw_pointer_cast(fill_.data()),
                      .prod       = thrust::raw_pointer_cast(prod_.data()),
                      .wmax       = thrust::raw_pointer_cast(wmax_.data()),
                      .emin       = thrust::raw_pointer_cast(emin_.data()),
                      .emax       = thrust::raw_pointer_cast(emax_.data()),
                      .nbin       = par_.nbin,
                      .ncell_edge = par_.ncell_edge,
                      .ncell      = ncell,
                      .nkeys      = nkeys_,
                      .keyed      = keyed_,
                      .mass       = mass_,
                      .llo        = llo_,
                      .dl         = dl_,
                      .mins       = mins_,
                      .maxs       = maxs_ };
}

void
  Pop::follow_container(const pic::ParticleContainer& con)
{
  // slots grown onto the container since the layout are tail, key -1
  const auto n = con.size();
  if(skey_.size() < n) { skey_.resize(n, std::int32_t { -1 }); }
}

Pop::Totals
  Pop::totals(const pic::ParticleContainer& con)
{
  // ---- one reduction over the slots of (count, weight, energy, momentum) of the live LPs
  const auto ids_mds = con.ids_mds();
  const auto vel_mds = con.vel_mds();
  const auto wgt_mds = con.wgt_mds();
  const auto weighted = con.has_particle_weights();
  const auto w0       = static_cast<double>(con.weight());
  const auto m        = con.args().mass;
  using V6            = toolbox::VecD<double, 6>;  // N, W, E, Px, Py, Pz of one slot
  const auto ord0     = thrust::counting_iterator<I>(0);
  // the six sums a live slot contributes, zero for a dead one
  const auto f        = [=](const I i) {
    if(ids_mds[i][] == runko::dead_prtc_id) { return V6 {}; }
    const double vx = vel_mds[i][0], vy = vel_mds[i][1], vz = vel_mds[i][2];
    const auto e    = sstd::sqrt(m * m + vx * vx + vy * vy + vz * vz);
    const auto wg   = weighted ? static_cast<double>(wgt_mds[i][]) : w0;
    return V6(1.0, wg, wg * e, wg * vx, wg * vy, wg * vz);
  };
  const auto w = tyvi::mdgrid_work {};
  // summed over the container
  const auto t = thrust::reduce(
    w.on_this(),
    thrust::make_transform_iterator(ord0, f),
    thrust::make_transform_iterator(ord0 + static_cast<std::ptrdiff_t>(con.size()), f),
    V6 {},
    [](const V6& a, const V6& b) { return a + b; });
  return Totals { t[0], t[1], t[2], { t[3], t[4], t[5] } };
}

std::vector<std::int32_t>
  Pop::get_slot_keys(const pic::ParticleContainer& con)
{
  this->follow_container(con);
  const auto n       = con.size();
  const auto ids_mds = con.ids_mds();
  // compact the keys of the alive slots on the device, then copy the short list to the host
  auto out           = thrust::device_vector<std::int32_t>(n);
  const auto w       = tyvi::mdgrid_work {};
  const auto end     = thrust::copy_if(
    w.on_this(),
    skey_.begin(),
    skey_.begin() + static_cast<std::ptrdiff_t>(n),
    thrust::counting_iterator<I>(0),
    out.begin(),
    [=](const I i) { return ids_mds[i][] != runko::dead_prtc_id; });
  const auto h = thrust::host_vector<std::int32_t>(out.begin(), end);
  return std::vector<std::int32_t>(h.begin(), h.end());
}

Pop::Tables
  Pop::get_tables()
{
  Tables t;
  // device tables to host copies, for tests and diagnostics
  thrust::host_vector<I> b = b_, cap = cap_, fill = fill_, prod = prod_;
  thrust::host_vector<float> wmax = wmax_, emin = emin_, emax = emax_;
  t.b.assign(b.begin(), b.end());
  t.cap.assign(cap.begin(), cap.end());
  t.fill.assign(fill.begin(), fill.end());
  t.prod.assign(prod.begin(), prod.end());
  t.wmax.assign(wmax.begin(), wmax.end());
  t.emin.assign(emin.begin(), emin.end());
  t.emax.assign(emax.begin(), emax.end());
  t.llo = llo_;
  t.dl  = dl_;
  return t;
}

}  // namespace qed
