# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""The QED diagnostics of a multi-rank run: every record, each rank sums its tiles' depth,
weights, energies, counts, budgets and the escaping spectra, and histograms its own photons
(weight per energy, weight distribution) and leptons (momentum); the sums are reduced to the
root, which writes qed_<lap>.npz and, at the end, qed_timeseries.npz."""

import os
import numpy as np
from mpi4py import MPI

SPECIES = ("ph", "e-", "e+")
X_BINS = np.logspace(-8, 2, 101)      # photon energy, for the spectra
P_BINS = np.logspace(-2, 3, 61)       # lepton momentum
W_BINS = np.logspace(-9, 1, 61)       # photon weight
MU_BINS = np.linspace(-1, 1, 33)


class QedDiagnostics:

    def __init__(self, simulation, conf, dt_box, n_tiles_edge, outdir):
        self.sim, self.conf, self.dt, self.nt, self.outdir = simulation, conf, dt_box, n_tiles_edge, outdir
        self.comm = MPI.COMM_WORLD
        self.series = {k: [] for k in ("lap", "t", "tau", "tau_min", "tau_max", "W", "E", "N", "ell_esc", "esc_W", "heat_E",
                                       "counts", "moved", "n_ev")}
        self.prev = None
        self.names = None

    def tile_sums(self):
        """per rank: [W x3, E x3, N x3, esc_E, esc_W, heat_E, inj_E, inj_W] and per-tile depths"""
        loc = np.zeros(14); taus = []
        counts = None; moved = None
        for tile in self.sim.local_tiles():
            for i, s in enumerate(SPECIES):
                tt = tile.totals(s); loc[i] += tt["W"]; loc[3 + i] += tt["E"]; loc[6 + i] += tt["N"]
            b = tile.budget()
            loc[9:14] += [b["esc_E"], b["esc_W"], b["heat_E"], b["inj_E"], b["inj_W"]]
            taus.append(tile.totals("e-")["W"] + tile.totals("e+")["W"])
            cn, mv = tile.counts(), tile.moved()
            if self.names is None:
                self.names = sorted(cn)
            c = np.array([cn[k] for k in self.names]); m = np.array([mv[k] for k in self.names])
            counts = c if counts is None else counts + c
            moved = m if moved is None else moved + m
        if counts is None:
            counts = np.zeros(len(self.names or [])); moved = np.zeros(len(self.names or []))
        return loc, np.array(taus), counts, moved

    def histograms(self):
        """per rank: escaping spectra (from the tiles, since the last record), photon and lepton histograms"""
        h = dict(esc=np.zeros(len(X_BINS) - 1), esc_mu=np.zeros((len(X_BINS) - 1, len(MU_BINS) - 1)),
                 ph_xw=np.zeros(len(X_BINS) - 1), ph_xn=np.zeros(len(X_BINS) - 1),
                 ph_wn=np.zeros(len(W_BINS) - 1), ph_xw2=np.zeros((len(X_BINS) - 1, len(W_BINS) - 1)),
                 lep_pw=np.zeros(len(P_BINS) - 1), lep_pn=np.zeros(len(P_BINS) - 1))
        for tile in self.sim.local_tiles():
            edges, hx = tile.escaped_spectrum()
            h["esc"] += hx; h["esc_mu"] += tile.escaped_mu_spectrum()
            tile.reset_budget()
            vx, vy, vz = (np.asarray(a) for a in tile.get_velocities(tile.species("ph")))
            x = np.sqrt(vx * vx + vy * vy + vz * vz); w = np.asarray(tile.get_weights(tile.species("ph")))
            h["ph_xw"] += np.histogram(x, X_BINS, weights=w * x)[0]
            h["ph_xn"] += np.histogram(x, X_BINS)[0]
            h["ph_wn"] += np.histogram(w, W_BINS)[0]
            h["ph_xw2"] += np.histogram2d(x, w, [X_BINS, W_BINS], weights=w)[0]
            for s in ("e-", "e+"):
                vx, vy, vz = (np.asarray(a) for a in tile.get_velocities(tile.species(s)))
                p = np.sqrt(vx * vx + vy * vy + vz * vz); w = np.asarray(tile.get_weights(tile.species(s)))
                h["lep_pw"] += np.histogram(p, P_BINS, weights=w)[0]
                h["lep_pn"] += np.histogram(p, P_BINS)[0]
        return h

    def record(self, lap, budget_span):
        """reduce this record to the root and write it; budget_span = laps since the last record.
        tau is the depth of the box, tau_min/max the depth the box would have at the least and
        the most loaded tile's pair density"""
        loc, taus, counts, moved = self.tile_sums()
        tot = self.comm.reduce(loc, op=MPI.SUM, root=0)
        tau_min = self.comm.reduce(taus.min() if len(taus) else np.inf, op=MPI.MIN, root=0)
        tau_max = self.comm.reduce(taus.max() if len(taus) else -np.inf, op=MPI.MAX, root=0)
        counts = self.comm.reduce(counts, op=MPI.SUM, root=0)
        moved = self.comm.reduce(moved, op=MPI.SUM, root=0)
        h = self.histograms()
        hs = {k: self.comm.reduce(v, op=MPI.SUM, root=0) for k, v in h.items()}
        if self.comm.rank != 0:
            return
        span = max(budget_span, 1) * self.dt
        nt2 = self.nt ** 2
        S = self.series
        S["lap"].append(lap); S["t"].append(lap * self.dt)
        S["tau"].append((tot[1] + tot[2]) / nt2); S["tau_min"].append(tau_min * self.nt); S["tau_max"].append(tau_max * self.nt)   # a tile's depth, scaled to the box edge
        S["W"].append(tot[0:3] / nt2); S["E"].append(tot[3:6] / nt2); S["N"].append(tot[6:9])
        S["ell_esc"].append(tot[9] / nt2 / span); S["esc_W"].append(tot[10] / nt2); S["heat_E"].append(tot[11] / nt2)
        dc = counts - (self.prev[0] if self.prev else 0); dm = moved - (self.prev[1] if self.prev else 0)
        self.prev = (counts, moved)
        S["counts"].append(dc / span); S["moved"].append(dm / nt2 / span); S["n_ev"].append(dc.sum())
        np.savez_compressed(os.path.join(self.outdir, f"qed_{lap}.npz"), lap=lap, t=lap * self.dt, span=span,
                            x_bins=X_BINS, p_bins=P_BINS, w_bins=W_BINS, mu_bins=MU_BINS,
                            esc=hs["esc"] / nt2 / span, esc_mu=hs["esc_mu"] / nt2 / span,
                            ph_xw=hs["ph_xw"] / nt2, ph_xn=hs["ph_xn"], ph_wn=hs["ph_wn"], ph_xw2=hs["ph_xw2"] / nt2,
                            lep_pw=hs["lep_pw"] / nt2, lep_pn=hs["lep_pn"], names=np.array(self.names))

    def save_series(self):
        if self.comm.rank != 0:
            return
        S = {k: np.array(v) for k, v in self.series.items()}
        S["names"] = np.array(self.names)
        np.savez_compressed(os.path.join(self.outdir, "qed_timeseries.npz"), **S)
