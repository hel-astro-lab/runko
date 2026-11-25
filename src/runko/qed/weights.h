// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// The weight bookkeeping of the QED Monte Carlo, as pure functions.
//
// Two LPs of weights w_i and w_j stand for w_i w_j 2 r / N_2 reactions per step, r
// being the rate function of the channel. An event moves the weight S and fires with
// the probability (w_i w_j / S) 2 r / N_2, so that probability times size is the
// reaction count whatever S is chosen. The incident is updated with the probability
// S / w_i; the partner takes S / w_j exchanges, a probability while below one and a
// multiplicity above.

#include "runko/tools/math.h"

namespace qed {

/// Exchanges a light partner of energy gam may take as one: kappa (gam - 1) / dgmax,
/// dgmax being the most energy a single exchange gives it; at least 1, and 1 for
/// kappa <= 0.
[[nodiscard]]
constexpr double
  fkin(const double kappa, const double gam, const double dgmax)
{
  if(kappa <= 0.0) { return 1.0; }
  return sstd::max(1.0, kappa * (gam - 1.0) / sstd::max(dgmax, 1e-300));
}

/// Weight an event of the pair moves: min(w_i, w_j) between equals (fk = 1), and up
/// to the whole of a heavy incident when the light partner takes fk exchanges as one.
[[nodiscard]]
constexpr double
  event_size(const double wi, const double wj, const double fk)
{ return sstd::min(wi, fk * wj); }

/// Weight of the pair in its event probability, w_i w_j / event_size. Rising in both
/// weights and falling in fk, so the largest weights of two groups and their smallest
/// fk bound it.
[[nodiscard]]
constexpr double
  pair_weight(const double wi, const double wj, const double fk)
{ return sstd::max(wj, wi / fk); }

}  // namespace qed
