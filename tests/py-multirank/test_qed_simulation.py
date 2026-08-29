# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""The QED tile on a 2x2x2 grid over several ranks: reactions at the top of the lap, photons
and their products pushed across tile boundaries and exchanged, the layout kept as tail
arrivals. Checks conservation across the exchange, unique ids, and every particle in its tile."""

import mpi_unittest
import tempfile
import numpy as np
from mpi4py import MPI

import runko

WL = 1e-3
NC = 8


def broadcast_outdir():
    comm = MPI.COMM_WORLD
    outdir = tempfile.mkdtemp(prefix="runko-qed-mpi-") if comm.rank == 0 else None
    return comm.bcast(outdir, root=0)


def make_grid(seed=1):
    config = runko.Configuration(None)
    config.tile_partitioning = "catepillar_track"
    config.catepillar_track_length = 1
    config.n_laps = 20
    config.n_tiles = [2, 2, 2]
    config.n_cells_per_tile = [NC, NC, NC]
    config.cfl = 0.45
    config.field_propagator = "fdtd2"
    config.q0, config.m0, config.w0 = -1, 1, WL
    config.q1, config.m1, config.w1 = 1, 1, WL
    config.q2, config.m2 = 0, 0
    config.particle_pusher = "boris"
    config.field_interpolator = "linear_1st"
    config.current_depositer = "zigzag_1st"
    config.qed_H_lamC = 4.111e9
    config.qed_nbin, config.qed_nsort, config.qed_slack, config.qed_smin = 8, 5, 0.2, 4
    config.qed_a_res = 0.0                      # no splitting: photon count and weight are conserved exactly
    config.qed_seed = seed
    config.qed_channels = ["compton-", "compton+"]
    config.io_outdir = broadcast_outdir()
    return config, runko.TileGrid(config)


def iso_dirs(n, rng):
    cz = -1 + 2 * rng.random(n); ph = 2 * np.pi * rng.random(n); sz = np.sqrt(1 - cz * cz)
    return np.stack([sz * np.cos(ph), sz * np.sin(ph), cz], 1)


def sample_mj(theta, n, rng):
    out = np.empty(n); k = 0
    while k < n:
        m = 2 * (n - k) + 10
        xi = rng.random((m, 4))
        eta = -theta * np.log(xi[:, 0] * xi[:, 1] * xi[:, 2]); zeta = -theta * np.log(np.prod(xi, 1))
        p = eta[zeta ** 2 - eta ** 2 > 1]
        take = min(len(p), n - k); out[k:k + take] = p[:take]; k += take
    return out


def rows(u):
    return [tuple(map(float, v)) for v in np.asarray(u).reshape(-1, 3)]


def totals(simulation):
    """summed over every tile of every rank: weight and count per species"""
    loc = np.zeros(6)
    for tile in simulation.local_tiles():
        for i, s in enumerate(("ph", "e-", "e+")):
            tt = tile.totals(s); loc[i] += tt["W"]; loc[3 + i] += tt["N"]
    return MPI.COMM_WORLD.allreduce(loc, op=MPI.SUM)


def test_qed_across_tiles():
    conf, tile_grid = make_grid()
    rank = MPI.COMM_WORLD.rank
    rng = np.random.default_rng(100 + rank)
    for idx in tile_grid.local_tile_indices():
        tile = runko.qed.threeD.Tile(idx, conf)
        zero = lambda x, y, z: (0, 0, 0)
        tile.set_EBJ(zero, zero, zero)
        nph, npm = 2000, 400
        x = rng.gamma(3.0, 0.05, nph)
        tile.set_particles("ph", rows(x[:, None] * iso_dirs(nph, rng)), list(WL * np.exp(rng.uniform(np.log(0.3), np.log(3.0), nph))))
        for s in ("e-", "e+"):
            p = sample_mj(2.0, npm, rng)
            tile.set_particles(s, rows(p[:, None] * iso_dirs(npm, rng)), [WL] * npm)
        tile_grid.add_tile(tile, idx)
    simulation = tile_grid.configure_simulation(conf)
    t0 = totals(simulation)

    def lap(x):
        x.prtcl_interact()
        x.prtcl_push()
        x.prtcl_pack_outgoing()
        x.comm_external(runko.tools.comm_mode.pic_particle)
        x.comm_local(runko.tools.comm_mode.pic_particle)

    simulation.for_each_lap(lap)
    t1 = totals(simulation)

    # weights and counts are conserved across all the exchanges (Compton moves no weight
    # between species and, without splitting, makes no new photon LPs)
    asserts = []
    for i in range(6):
        asserts.append(mpi_unittest.assertEqualDeferred(True, bool(abs(t1[i] - t0[i]) <= 1e-9 * max(t0[i], 1.0))))
    # events happened, on every rank
    n_ev = sum(sum(tile.counts().values()) for tile in simulation.local_tiles())
    asserts.append(mpi_unittest.assertEqualDeferred(True, n_ev > 0))
    # every particle inside its own tile, with the layout intact
    for tile in simulation.local_tiles():
        mins, maxs = np.array(tile.mins), np.array(tile.maxs)
        for s in (0, 1, 2):
            pos = np.stack([np.asarray(a) for a in tile.get_positions(s)], 1)
            inside = np.all((pos >= mins) & (pos < maxs)) if len(pos) else True
            asserts.append(mpi_unittest.assertEqualDeferred(True, bool(inside)))
        asserts.append(mpi_unittest.assertEqualDeferred(True, tile.laid()))
    mpi_unittest.assertDeferredResults(asserts)
    # ids unique per species over the whole grid
    for s in (0, 1, 2):
        ids = np.concatenate([np.asarray(tile.get_ids(s)) for tile in simulation.local_tiles()] or [np.zeros(0)])
        all_ids = MPI.COMM_WORLD.gather(ids, root=0)
        if rank == 0:
            a = np.concatenate(all_ids)
            mpi_unittest.assertEqual(len(a), len(np.unique(a)))
        else:
            mpi_unittest.assertEqual(0, 0)


if __name__ == "__main__":
    test_qed_across_tiles()
